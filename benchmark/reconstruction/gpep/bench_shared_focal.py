#!/usr/bin/env python3
# Copyright (c), ETH Zurich and UNC Chapel Hill.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#
#     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
#       its contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""End-to-end ETH3D benchmark for the shared-focal two-view geometry path.

Reconstructs the same scenes twice -- once with a baseline colmap binary and
once with the branch under test -- and prints the paired, per-seed comparison
of the two. Both arms share one set of scene workspaces, so features and raw
matches are extracted once and reused.

WHY THIS SCRIPT EXISTS

  Nothing about a stock benchmark run reaches the shared-focal estimator.
  `EstimateTwoViewGeometry` only dispatches to it when the two images share a
  camera and that camera has no focal prior, so the benchmark has to be run
  with both of these:

    --uncalibrated   stops the harness from priming every camera with its GT
                     focal, which would select the calibrated (essential
                     matrix) path instead.
    --single_camera  makes every image of a scene share one camera. Without
                     it the extractor gives ETH3D one camera per image and no
                     pair ever shares one.

  The harness has no passthrough for `automatic_reconstructor` options, so
  this script wraps each colmap binary in a small generated shell shim.
  The shim appends `--single_camera 1`, and just before the matching
  phase it also clears the focal prior outright and restores the
  extractor's default focal. That last step is belt and braces:
  --uncalibrated only means "do not write GT intrinsics", so a prior left
  in a reused workspace by an earlier run would otherwise survive and
  route every pair to the calibrated path with nothing to show for it.
  Nothing under benchmark/ is modified; the shims are written into the run
  directory and are readable there.

SCENE SELECTION

  Only the ETH3D DSLR scenes whose ground truth has exactly one camera are
  evaluated by default -- 19 of 25, everything except courtyard, electro,
  facade, playground, terrace and terrains. Forcing a single camera on a
  scene that was genuinely shot with several models a calibration the ground
  truth contradicts, so those scenes would measure the wrong thing. The list
  is recomputed from the ground truth on every run rather than hard-coded.

  This restriction is specific to forcing a shared camera. Its sibling
  bench_one_sided_focal.py asserts nothing about the physical camera and runs
  on all 25 scenes, so the two benchmarks are not comparable scene for scene.

READING THE RESULT

  ETH3D is bistable on several scenes: a single seed says very little, and a
  difference smaller than the printed standard deviation is noise. The
  original evaluation of the 6pt shared-focal work needed 5 seeds per arm to
  resolve a +0.9 AUC@1 effect, and per-scene deltas from one seed were
  entirely misleading. Prefer more seeds over more scenes.

PREREQUISITES

  * the ETH3D DSLR data: `python ../download.py --datasets eth3d`
  * pycolmap importable: `pip install -r ../requirements.txt`
  * two colmap binaries, ideally built from the same source tree at two
    different commits, with identical cmake options

EXAMPLES

  # Smoke test on one small scene, one seed (~minutes):
  python bench_shared_focal.py --baseline-colmap A --branch-colmap B \\
      --scenes door --num-seeds 1 --run-name shared-focal-smoke

  # Real run (hours):
  python bench_shared_focal.py --baseline-colmap A --branch-colmap B \\
      --num-seeds 5

  # Re-print the comparison of a finished run:
  python bench_shared_focal.py --baseline-colmap A --branch-colmap B \\
      --compare-only
"""

import argparse
import shlex
import sqlite3
import stat
import subprocess
import sys
from pathlib import Path

import pycolmap

# These drivers live one level below the harness they call into.
HARNESS_DIR = Path(__file__).resolve().parent.parent

# The shared-focal setting is only defined for the DSLR category; the ETH3D
# rig scenes are multi-sensor by construction.
CATEGORY = "dslr"

# Logged by colmap when the SIFT GPU matcher cannot allocate. The scene
# then keeps zero matches without any other sign of trouble.
GPU_OOM_MARKER = "Not enough GPU memory"

# Focal length the extractor assigns to an image with no EXIF focal, as a
# multiple of its longer side. See ImageReaderOptions.
DEFAULT_FOCAL_LENGTH_FACTOR = 1.2

SHIM_TEMPLATE = r"""#!/usr/bin/env bash
# Generated by bench_shared_focal.py. Rewritten on every run; do not edit.
#
# Two things the harness cannot express, both of them on
# `colmap automatic_reconstructor`, which it offers no passthrough for:
#
#   * --single_camera 1 makes every image of a scene share one camera.
#   * Just before the matching phase, the shared camera's focal prior is
#     cleared and its focal reset to the extractor's default. --uncalibrated
#     only stops the harness from writing GT intrinsics; it cannot undo a
#     prior left behind in a workspace by an earlier, differently configured
#     run, and such a prior would silently disable the shared-focal path.
#
# Together these select the shared-focal two-view geometry path. Other
# subcommands (database_cleaner) are forwarded untouched.
set -euo pipefail
if [ "${1-}" = automatic_reconstructor ]; then
  matching=0 workspace= previous= phase=other seed=-1
  for argument in "$@"; do
    case $previous in
      --matching) matching=$argument ;;
      --random_seed) seed=$argument ;;
      --workspace_path) workspace=$argument ;;
    esac
    previous=$argument
  done
  # The harness runs one phase per invocation, flagged 1; name it for the
  # timing row.
  for candidate in extraction matching sparse; do
    case " $* " in *" --$candidate 1 "*) phase=$candidate ;; esac
  done
  if [ "$matching" = 1 ]; then
    @PYTHON@ @DRIVER@ --clear-focal-priors \
      --database "$workspace/database.db"
  fi
  # Time the phase rather than exec'ing into it. colmap logs its own elapsed
  # time, but the harness overwrites each scene's phase log on every seed, so
  # only the final seed of the final arm would survive to be compared.
  started=$(date +%s)
  set +e
  @COLMAP@ "$@" --single_camera 1
  status=$?
  set -e
  elapsed=$(($(date +%s) - started))
  printf '%s,%s,%s,%s\n' "$(basename "$workspace")" "$phase" "$seed" \
    "$elapsed" >> @TIMING@
  exit $status
fi
exec @COLMAP@ "$@"
"""


def single_camera_scenes(data_path: Path) -> list[str]:
    """Returns the ETH3D DSLR scenes whose GT model has exactly one camera."""
    root = data_path / "eth3d" / CATEGORY
    if not root.is_dir():
        raise SystemExit(
            f"No ETH3D {CATEGORY} data under {root}.\n"
            f"Get it with: python {HARNESS_DIR / 'download.py'} "
            "--datasets eth3d"
        )
    scenes = []
    for scene_path in sorted(p for p in root.iterdir() if p.is_dir()):
        cameras_txt = next(
            iter(sorted(scene_path.glob("*_calibration*/cameras.txt"))), None
        )
        if cameras_txt is None:
            continue
        num_cameras = sum(
            1
            for line in cameras_txt.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        )
        if num_cameras == 1:
            scenes.append(scene_path.name)
    if not scenes:
        raise SystemExit(f"No single-camera scenes found under {root}")
    return scenes


def clear_focal_priors(database_path: Path) -> None:
    """Drops every focal prior in the database and restores default focals.

    Written explicitly rather than relying on --uncalibrated, which only means
    "do not write GT intrinsics": workspaces are shared and long-lived, so a
    prior from an earlier run would otherwise survive and route every pair to
    the calibrated path instead.
    """
    with pycolmap.Database.open(str(database_path)) as database:
        cameras = database.read_all_cameras()
        for camera in cameras:
            camera.focal_length = DEFAULT_FOCAL_LENGTH_FACTOR * max(
                camera.width, camera.height
            )
            camera.has_prior_focal_length = False
            database.update_camera(camera)
    print(
        f"[bench_shared_focal] {database_path.parent.name}: cleared the "
        f"focal prior of {len(cameras)} camera(s)",
        flush=True,
    )


def write_shim(shim_path: Path, colmap_path: Path) -> Path:
    """Writes the executable wrapper around one colmap binary."""
    shim_path.parent.mkdir(parents=True, exist_ok=True)
    shim = SHIM_TEMPLATE
    for placeholder, value in (
        ("@PYTHON@", sys.executable),
        ("@DRIVER@", str(Path(__file__).resolve())),
        ("@COLMAP@", str(colmap_path)),
        ("@TIMING@",
         str(timing_path(shim_path.parent.parent, shim_path.stem))),
    ):
        shim = shim.replace(placeholder, shlex.quote(value))
    shim_path.write_text(shim)
    shim_path.chmod(shim_path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP)
    # The shim appends one row per colmap call, so start each arm from an
    # empty file: otherwise rows from an earlier or interrupted run of the
    # same arm would be folded into its runtime report.
    timing_path(shim_path.parent.parent, shim_path.stem).write_text("")
    return shim_path


def run_arm(
    args: argparse.Namespace, label: str, colmap_path: Path, scenes: list[str]
) -> None:
    """Evaluates one binary over every seed, writing <label>_s<seed>.pkl."""
    run_dir = args.run_path / args.run_name
    shim = write_shim(run_dir / ".shims" / f"{label}.sh", colmap_path)
    command = [
        sys.executable,
        str(HARNESS_DIR / "evaluate.py"),
        "--colmap_path", str(shim),
        "--report_name", label,
        "--num_seeds", str(args.num_seeds),
        # RANSAC seeds per thread as random_seed + omp_get_thread_num(), so
        # only single-threaded scenes are reproducible at a fixed seed. Use
        # --num-parallel-scenes for throughput instead: parallelism across
        # scenes does not touch any one scene's random number stream.
        "--threads_per_scene", str(args.threads_per_scene),
        "--num_parallel_scenes", str(args.num_parallel_scenes),
        "--gpu_index", args.gpu_index,
        "--run_path", str(args.run_path),
        "--run_name", args.run_name,
        "--data_path", str(args.data_path),
        "--datasets", "eth3d",
        "--categories", CATEGORY,
        "--scenes", *scenes,
        "--mapper", args.mapper,
        "--feature", "sift",
        "--quality", args.quality,
        # Necessary but not sufficient: it only stops the harness from
        # writing GT intrinsics. The shim clears the prior outright.
        "--uncalibrated",
        f"--overwrite_{args.overwrite}",
        "--progress" if args.progress else "--no-progress",
        *args.extra_eval_arg,
    ]
    print(f"\n==== arm '{label}': {colmap_path} ====", flush=True)
    print(shlex.join(command), flush=True)
    if subprocess.run(command, check=False).returncode != 0:
        raise SystemExit(
            f"Arm '{label}' failed. The harness swallows colmap's "
            f"output into per-phase logs, so look at "
            f"{run_dir / 'eth3d' / CATEGORY}/<scene>/"
            "{extraction,matching,reconstruction}.log -- which is also "
            "where this script's own per-scene line is written."
        )


def summarize_databases(args: argparse.Namespace, scenes: list[str]) -> None:
    """Prints what actually landed in each scene database.

    Cheap insurance against a run that looks like a regression but is really a
    misconfiguration: every scene must show exactly 1 camera and 0 cameras
    with a focal prior, or the shared-focal path was never taken. A scene with
    no verified pairs matched nothing at all, which on a shared workspace is
    almost always GPU exhaustion from too many parallel scenes.
    """
    run_dir = args.run_path / args.run_name / "eth3d" / CATEGORY
    print(f"\n==== database summary: {run_dir} ====")
    print(f"{'scene':<20}{'images':>8}{'cameras':>9}{'priors':>8}{'pairs':>9}")
    empty_scenes = []
    for scene in scenes:
        database_path = run_dir / scene / "database.db"
        if not database_path.exists():
            print(f"{scene:<20}{'-- no database --':>34}")
            continue
        with sqlite3.connect(f"file:{database_path}?mode=ro", uri=True) as db:
            count = lambda query: db.execute(query).fetchone()[0]  # noqa: E731
            num_images = count("SELECT COUNT(*) FROM images")
            num_cameras = count("SELECT COUNT(*) FROM cameras")
            num_priors = count(
                "SELECT COUNT(*) FROM cameras WHERE prior_focal_length != 0"
            )
            num_pairs = count(
                "SELECT COUNT(*) FROM two_view_geometries WHERE rows > 0"
            )
        if num_pairs == 0:
            empty_scenes.append(scene)
        flag = "" if num_cameras == 1 and num_priors == 0 else "  <-- CHECK"
        print(
            f"{scene:<20}{num_images:>8}{num_cameras:>9}"
            f"{num_priors:>8}{num_pairs:>9}{flag}"
        )
    print(
        "\nExpected: 1 camera and 0 priors per scene. Anything else means the "
        "shared-focal\npath was not exercised -- most likely a workspace "
        "reused from a differently\nconfigured run, which needs a fresh "
        "--run-name or --overwrite database."
    )
    abort_on_empty_scenes(run_dir, empty_scenes)


def abort_on_empty_scenes(run_dir: Path, empty_scenes: list[str]) -> None:
    """Aborts when a scene ended up with no verified pairs.

    Nearly always GPU exhaustion from too many parallel scenes: colmap logs
    "Not enough GPU memory to match N features", the matcher worker never
    starts, and geometric verification writes nothing. The raw matches survive
    -- it is the two-view geometries, just cleared by --overwrite, that are
    never regenerated -- so the scene registers no images and its AUC
    collapses. The collapse lands on whichever arm happened to be running, so
    it does not cancel out of the paired difference and reads as a large
    regression. Failing here is the difference between rerunning and quoting a
    wrong number.
    """
    if not empty_scenes:
        return
    starved = [
        scene
        for scene in empty_scenes
        if (log := run_dir / scene / "matching.log").exists()
        and GPU_OOM_MARKER in log.read_text(errors="replace")
    ]
    detail = (
        f"{len(starved)} of them logged {GPU_OOM_MARKER!r}: "
        f"{' '.join(starved)}"
        if starved
        else "None of them logged a GPU memory error, so check their "
        "matching.log for another cause"
    )
    raise SystemExit(
        f"\nABORT: {len(empty_scenes)} scene(s) have no verified pairs and "
        f"reconstructed nothing:\n  {' '.join(empty_scenes)}\n{detail}.\n"
        "Their numbers are meaningless and the comparison is not printed. "
        "Lower\n--num-parallel-scenes (1 is always safe) and run the same "
        "command again;\naffected scenes are re-matched automatically."
    )


def timing_path(run_dir: Path, label: str) -> Path:
    """Where the shim appends one `scene,phase,seconds` row per colmap call."""
    return run_dir / f"timing-{label}.csv"


def report_runtime(run_dir: Path, labels: tuple[str, str]) -> None:
    """Prints the cost of geometric verification, normalised per image pair.

    Seconds per scene are not comparable -- scenes differ by two orders of
    magnitude in pair count -- so the figure reported is milliseconds per
    pair that verification actually ran on, which is every pair holding raw
    matches. Totals are summed before dividing, so large scenes weigh
    proportionally rather than each scene counting equally.

    This is the phase total, not solver time. At --quality high the phase also
    re-matches every pair under the estimated geometry, and that dominates, so
    a slowdown confined to the solver is diluted here. It is nonetheless what
    the pipeline actually pays. Isolating the estimator wants a microbenchmark.

    Results are broken down by seed because the very first arm-seed run
    against a cold workspace also pays for raw matching, which inflates it;
    that row is visible rather than silently averaged in.
    """
    pairs_by_scene = {}
    category_dir = run_dir / "eth3d" / CATEGORY
    if category_dir.is_dir():
        scenes = (p for p in category_dir.iterdir() if p.is_dir())
        for scene_path in sorted(scenes):
            database_path = scene_path / "database.db"
            if not database_path.exists():
                continue
            with sqlite3.connect(
                f"file:{database_path}?mode=ro", uri=True
            ) as database:
                pairs_by_scene[scene_path.name] = database.execute(
                    "SELECT COUNT(*) FROM matches WHERE rows > 0"
                ).fetchone()[0]

    runs: dict[str, dict[str, tuple[int, int]]] = {}
    for label in labels:
        path = timing_path(run_dir, label)
        if not path.exists():
            print(f"\nNo runtime recorded for arm '{label}'")
            return
        by_seed: dict[str, tuple[int, int]] = {}
        for line in path.read_text().splitlines():
            fields = line.split(",")
            if len(fields) != 4 or fields[1] != "matching":
                continue
            scene, _, seed, seconds = fields
            pairs = pairs_by_scene.get(scene)
            if not pairs:
                continue
            total_seconds, total_pairs = by_seed.get(seed, (0, 0))
            by_seed[seed] = (total_seconds + int(seconds), total_pairs + pairs)
        if not by_seed:
            print(f"\nNo runtime recorded for arm '{label}'")
            return
        runs[label] = by_seed

    print("\n==== geometric verification cost (matching phase) ====")
    print(f"{'arm':<14}{'seed':>6}{'pairs':>10}{'total':>10}{'ms/pair':>10}")
    overall = {}
    for label in labels:
        seconds_sum = pairs_sum = 0
        for seed in sorted(runs[label], key=int):
            seconds, pairs = runs[label][seed]
            seconds_sum += seconds
            pairs_sum += pairs
            print(
                f"{label:<14}{seed:>6}{pairs:>10}"
                f"{seconds / 60:>9.1f}m{1000 * seconds / pairs:>9.1f}"
            )
        overall[label] = 1000 * seconds_sum / pairs_sum
    branch, baseline = labels
    change = 100.0 * (overall[branch] - overall[baseline]) / overall[baseline]
    print(
        f"\n{branch}: {overall[branch]:.1f} ms/pair, "
        f"{baseline}: {overall[baseline]:.1f} ms/pair  ->  {change:+.1f}%"
    )


def compare_arms(args: argparse.Namespace) -> None:
    """Prints the branch, baseline and paired branch-baseline tables.

    Only the seeds present in both arms are compared.
    """
    run_dir = args.run_path / args.run_name
    # compare.py prints A - B, so the branch goes in as A: a positive delta
    # then means the branch is ahead of the baseline, which is the direction
    # anyone reading an A/B expects.
    command = [
        sys.executable,
        str(HARNESS_DIR / "compare.py"),
        "--report_a_path_prefix", str(run_dir / args.branch_label),
        "--report_b_path_prefix", str(run_dir / args.baseline_label),
        "--labels", args.branch_label, args.baseline_label,
    ]
    print(f"\n==== comparison ====\n{shlex.join(command)}", flush=True)
    subprocess.run(command, check=True)


def parse_clear_args(argv: list[str]) -> argparse.Namespace:
    """Parses the shim-only invocation that clears one scene's priors."""
    parser = argparse.ArgumentParser(
        prog=f"{Path(__file__).name} --clear-focal-priors",
        description=clear_focal_priors.__doc__,
    )
    parser.add_argument("--clear-focal-priors", action="store_true")
    parser.add_argument("--database", type=Path, required=True)
    return parser.parse_args(argv)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--baseline-colmap",
        type=Path,
        required=True,
        metavar="BIN",
        help="Colmap binary for the A arm, usually built from main.",
    )
    parser.add_argument(
        "--branch-colmap",
        type=Path,
        required=True,
        metavar="BIN",
        help="Colmap binary for the B arm, built from the branch under test.",
    )
    parser.add_argument("--baseline-label", default="baseline")
    parser.add_argument("--branch-label", default="branch")
    parser.add_argument(
        "--data-path",
        type=Path,
        default=HARNESS_DIR / "data",
        help="Dataset root that download.py populated.",
    )
    parser.add_argument(
        "--run-path",
        type=Path,
        default=HARNESS_DIR / "runs",
        help="Root holding the scene workspaces and reports of all runs.",
    )
    parser.add_argument(
        "--run-name",
        default="shared-focal",
        help="Subdirectory of --run-path shared by both arms. Use a fresh "
        "name whenever the configuration changes, since a workspace caches "
        "the database it was built with.",
    )
    parser.add_argument(
        "--scenes",
        nargs="+",
        default=None,
        metavar="SCENE",
        help="Scenes to evaluate. Defaults to every ETH3D DSLR scene whose "
        "GT has a single camera.",
    )
    parser.add_argument(
        "--num-seeds",
        type=int,
        default=3,
        help="Random seeds per arm, run as seeds 0..N-1. The seeds are "
        "paired between the arms, so A-B is a paired difference.",
    )
    parser.add_argument(
        "--num-parallel-scenes",
        type=int,
        default=1,
        help="Scenes reconstructed concurrently, and the only throughput "
        "knob that does not affect the result. Each concurrent scene runs "
        "two SIFT GPU matchers, and exhausting GPU memory makes the "
        "affected scenes keep zero matches and register nothing -- which "
        "reads like an algorithmic regression but is not. On ETH3D DSLR "
        "at --quality high, 3 was already enough to starve a 10 GB card "
        "intermittently, so raise this only while watching for the abort "
        "at the end of a run. Lower it rather than switching to the CPU.",
    )
    parser.add_argument(
        "--threads-per-scene",
        type=int,
        default=1,
        help="Threads within one scene. Leave at 1: RANSAC seeds per thread, "
        "so a fixed seed is only reproducible single-threaded, and the "
        "paired comparison relies on that.",
    )
    parser.add_argument("--gpu-index", default="-1")
    parser.add_argument(
        "--mapper",
        default="incremental",
        choices=["incremental", "hierarchical", "global"],
    )
    parser.add_argument(
        "--quality", default="high", choices=["low", "medium", "high"]
    )
    parser.add_argument(
        "--overwrite",
        default="two_view_geometries",
        choices=["two_view_geometries", "matches", "database"],
        help="What each arm recomputes in the shared workspaces. "
        "two_view_geometries is both the cheapest and the fairest: raw "
        "matches are identical for the two arms and are not rewritten by "
        "geometric verification, so both arms verify exactly the same input.",
    )
    parser.add_argument(
        "--extra-eval-arg",
        action="append",
        default=[],
        metavar="ARG",
        help="Extra argument forwarded verbatim to evaluate.py. Repeatable.",
    )
    parser.add_argument(
        "--progress",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Live display of in-flight scenes. Off by default, since these "
        "runs are usually logged to a file.",
    )
    parser.add_argument(
        "--compare-only",
        action="store_true",
        help="Skip both arms and only re-print the comparison of a run.",
    )
    args = parser.parse_args()
    # The harness runs colmap with cwd set to the scene workspace, so every
    # path it is given has to be absolute.
    args.data_path = args.data_path.resolve()
    args.run_path = args.run_path.resolve()
    if args.baseline_label == args.branch_label:
        parser.error("--baseline-label and --branch-label must differ")
    return args


def main() -> None:
    if "--clear-focal-priors" in sys.argv[1:]:
        clear_focal_priors(parse_clear_args(sys.argv[1:]).database)
        return

    args = parse_args()
    scenes = args.scenes or single_camera_scenes(args.data_path)

    if args.compare_only:
        compare_arms(args)
        return

    for path in (args.baseline_colmap, args.branch_colmap):
        if not path.is_file():
            raise SystemExit(f"Not a colmap binary: {path}")
    if args.baseline_colmap.resolve() == args.branch_colmap.resolve():
        print("WARNING: both arms use the same binary; A-B measures only "
              "run-to-run noise.")

    print(f"Scenes ({len(scenes)}): {' '.join(scenes)}")
    print(f"Seeds: 0..{args.num_seeds - 1}")
    print(f"Run directory: {args.run_path / args.run_name}")

    # Sequentially: the arms share the scene workspaces and would corrupt each
    # other's databases if run at the same time.
    run_arm(args, args.baseline_label, args.baseline_colmap, scenes)
    run_arm(args, args.branch_label, args.branch_colmap, scenes)

    summarize_databases(args, scenes)
    report_runtime(
        args.run_path / args.run_name,
        (args.branch_label, args.baseline_label),
    )
    compare_arms(args)


if __name__ == "__main__":
    main()
