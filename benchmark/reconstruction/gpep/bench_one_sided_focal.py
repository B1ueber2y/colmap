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

"""End-to-end ETH3D benchmark for the one-sided focal two-view geometry path.

Reconstructs the same scenes twice -- once with a baseline colmap binary and
once with the branch under test -- and prints the paired, per-seed comparison
of the two. Both arms share one set of scene workspaces, so features and raw
matches are extracted once and reused.

WHY THIS SCRIPT EXISTS

  No stock benchmark configuration reaches the one-sided focal estimator.
  `EstimateTwoViewGeometry` dispatches to it only for a pair with exactly one
  calibrated side, and the harness primes either every camera with its GT
  focal or none at all: ETH3D calibrated gives 100% calibrated pairs, ETH3D
  --uncalibrated and IMC give 100% uncalibrated pairs, and zero mixed pairs
  either way. An A/B run without the setup below therefore proves nothing,
  however good it looks.

  So the harness is run with --uncalibrated -- it then never touches the
  cameras -- and each colmap binary is wrapped in a generated shell shim that
  gives a deterministic fraction of each scene's cameras their ground-truth
  focal right before the matching phase. With the default fraction of 0.5,
  about half of all pairs have exactly one calibrated side. Nothing under
  benchmark/ is modified; the shims are written into the run directory and
  are readable there.

  This works on ETH3D because its DSLR images carry no EXIF focal, so the
  extractor gives every image its own camera and calibrating a subset needs
  no camera surgery. Those cameras get the full GT intrinsics, exactly as the
  harness's own calibrated mode writes them.

SCENE SELECTION

  All 25 ETH3D DSLR scenes, and no selection is needed. One camera per image
  is not a property of a scene but of the extractor, which behaves that way
  on every one of them, so there is nothing to exclude. The split is per
  image as well, which is why a scene with several GT cameras is fine: each
  calibrated image is given the intrinsics of its own GT camera, and nothing
  is asserted that the ground truth contradicts.

  Its sibling bench_shared_focal.py does restrict to 19 scenes, because
  --single_camera asserts a single physical camera per scene, which is only
  true for some of them.

READING THE RESULT

  Expect e2e to be a no-regression check rather than the headline number.
  Bundle adjustment self-calibrates, so post-BA focal error lands around
  0.2% whether or not the pair initialization recovered the focal; the large
  effect of this path (two-view pose AUC@5 and focal error) is visible at the
  pair level and mostly erased by the time a reconstruction is finished.

  ETH3D is also bistable on several scenes: a single seed says very little,
  and a difference smaller than the printed standard deviation is noise.
  Prefer more seeds over more scenes.

PREREQUISITES

  * the ETH3D DSLR data: `python ../download.py --datasets eth3d`
  * pycolmap importable: `pip install -r ../requirements.txt`
  * two colmap binaries, ideally built from the same source tree at two
    different commits, with identical cmake options

EXAMPLES

  # Smoke test on one small scene, one seed (~minutes):
  python bench_one_sided_focal.py --baseline-colmap A --branch-colmap B \\
      --scenes door --num-seeds 1 --run-name one-sided-smoke

  # Real run (hours):
  python bench_one_sided_focal.py --baseline-colmap A --branch-colmap B \\
      --num-seeds 5

  # Re-print the comparison of a finished run:
  python bench_one_sided_focal.py --baseline-colmap A --branch-colmap B \\
      --compare-only
"""

import argparse
import random
import shlex
import sqlite3
import stat
import subprocess
import sys
from pathlib import Path

import pycolmap

# These drivers live one level below the harness they call into.
HARNESS_DIR = Path(__file__).resolve().parent.parent

CATEGORY = "dslr"

# Logged by colmap when the SIFT GPU matcher cannot allocate. The scene
# then keeps zero matches without any other sign of trouble.
GPU_OOM_MARKER = "Not enough GPU memory"

# Focal length the extractor assigns to an image with no EXIF focal, as a
# multiple of its longer side. See ImageReaderOptions.
DEFAULT_FOCAL_LENGTH_FACTOR = 1.2

# Image pair ids are image_id1 * kMaxNumImages + image_id2, see
# src/colmap/util/types.h.
K_MAX_NUM_IMAGES = 2147483647

SHIM_TEMPLATE = r"""#!/usr/bin/env bash
# Generated by bench_one_sided_focal.py. Rewritten on every run; do not edit.
#
# The harness primes either all of a scene's cameras with their GT focal or
# none of them, so it never produces a pair with exactly one calibrated side.
# It is run here with --uncalibrated, which primes nothing, and this wrapper
# writes the intended split itself just before the matching phase -- the point
# where the database exists and geometric verification has not run yet. The
# harness offers no passthrough for `automatic_reconstructor`, hence a
# wrapper. Other subcommands (database_cleaner) are forwarded untouched.
set -euo pipefail
if [ "${1-}" = automatic_reconstructor ]; then
  matching=0 workspace= images= previous= phase=other seed=-1
  for argument in "$@"; do
    case $previous in
      --matching) matching=$argument ;;
      --random_seed) seed=$argument ;;
      --workspace_path) workspace=$argument ;;
      --image_path) images=$argument ;;
    esac
    previous=$argument
  done
  # The harness runs one phase per invocation, flagged 1; name it for the
  # timing row.
  for candidate in extraction matching sparse; do
    case " $* " in *" --$candidate 1 "*) phase=$candidate ;; esac
  done
  if [ "$matching" = 1 ]; then
    @PYTHON@ @DRIVER@ --set-calibration-split \
      --database "$workspace/database.db" \
      --image-path "$images" \
      --calibrated-fraction @FRACTION@ \
      --split-seed @SPLIT_SEED@
  fi
  # Time the phase rather than exec'ing into it. colmap logs its own elapsed
  # time, but the harness overwrites each scene's phase log on every seed, so
  # only the final seed of the final arm would survive to be compared.
  started=$(date +%s)
  set +e
  @COLMAP@ "$@"
  status=$?
  set -e
  elapsed=$(($(date +%s) - started))
  printf '%s,%s,%s,%s\n' "$(basename "$workspace")" "$phase" "$seed" \
    "$elapsed" >> @TIMING@
  exit $status
fi
exec @COLMAP@ "$@"
"""


def ground_truth_path(image_path: Path) -> Path:
    """Locates an ETH3D scene's GT model from its image directory."""
    candidates = sorted(image_path.parent.glob("*_calibration*"))
    if not candidates:
        raise SystemExit(f"No GT calibration next to {image_path}")
    return candidates[0]


def set_calibration_split(
    database_path: Path, image_path: Path, fraction: float, seed: int
) -> None:
    """Splits a scene's cameras into a calibrated and an uncalibrated half.

    The calibrated ones are given the full GT intrinsics and a focal prior,
    exactly as the harness's own calibrated mode writes them. The rest are
    reset to the extractor's default focal with no prior, which is what a
    genuinely unknown camera looks like.

    Both halves are written explicitly rather than left as extraction produced
    them: workspaces are shared and long-lived, and one that an earlier,
    differently configured run already primed would otherwise skew the split
    without any visible sign.

    Selection is over cameras rather than images, so a camera shared by
    several images can never end up half primed, and it is drawn from a fixed
    generator so that both arms and every seed see the same partition.
    """
    ground_truth = pycolmap.Reconstruction(str(ground_truth_path(image_path)))
    gt_images_by_name = {
        image.name: image for image in ground_truth.images.values()
    }

    with pycolmap.Database.open(str(database_path)) as database:
        image_by_camera_id: dict[int, pycolmap.Image] = {}
        for image in sorted(database.read_all_images(), key=lambda i: i.name):
            image_by_camera_id.setdefault(image.camera_id, image)
        camera_ids = sorted(image_by_camera_id)

        # Shuffled rather than sliced: neighbouring images are the covisible
        # ones, so calibrating a contiguous half would confine the mixed pairs
        # to the seam between the two halves.
        num_calibrated = round(len(camera_ids) * fraction)
        calibrated = set(
            random.Random(seed).sample(camera_ids, num_calibrated)
        )

        missing = 0
        for camera in database.read_all_cameras():
            gt_image = None
            if camera.camera_id in calibrated:
                image = image_by_camera_id.get(camera.camera_id)
                gt_image = gt_images_by_name.get(image.name) if image else None
                if gt_image is None:
                    missing += 1
            if gt_image is None:
                camera.focal_length = DEFAULT_FOCAL_LENGTH_FACTOR * max(
                    camera.width, camera.height
                )
                camera.has_prior_focal_length = False
                database.update_camera(camera)
            else:
                gt_camera = ground_truth.cameras[gt_image.camera_id]
                gt_camera.camera_id = camera.camera_id
                gt_camera.has_prior_focal_length = True
                database.update_camera(gt_camera)

    print(
        f"[bench_one_sided_focal] {database_path.parent.name}: calibrated "
        f"{len(calibrated) - missing} of {len(camera_ids)} cameras"
        + (f" ({missing} absent from GT)" if missing else ""),
        flush=True,
    )


def write_shim(
    shim_path: Path, colmap_path: Path, args: argparse.Namespace
) -> Path:
    """Writes the executable wrapper that applies the calibration split."""
    shim_path.parent.mkdir(parents=True, exist_ok=True)
    shim = SHIM_TEMPLATE
    for placeholder, value in (
        ("@PYTHON@", sys.executable),
        ("@DRIVER@", str(Path(__file__).resolve())),
        ("@FRACTION@", str(args.calibrated_fraction)),
        ("@SPLIT_SEED@", str(args.split_seed)),
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
    shim = write_shim(run_dir / ".shims" / f"{label}.sh", colmap_path, args)
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
        *(["--scenes", *scenes] if scenes else []),
        "--mapper", args.mapper,
        "--feature", "sift",
        "--quality", args.quality,
        # The shim writes the split; the harness must not prime anything.
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


def summarize_databases(args: argparse.Namespace) -> None:
    """Prints the calibration mix each scene database ended up with.

    This is the check that the run measured anything at all: 'mixed'
    counts the verified pairs with exactly one calibrated side, and those
    are the only pairs the one-sided focal estimator ever sees. If it is
    zero everywhere, the two arms cannot differ whatever the branch does.
    """
    run_dir = args.run_path / args.run_name / "eth3d" / CATEGORY
    if not run_dir.is_dir():
        print(f"\nNo workspaces under {run_dir}")
        return
    print(f"\n==== database summary: {run_dir} ====")
    print(f"{'scene':<20}{'cameras':>9}{'calib':>8}"
          f"{'pairs':>9}{'mixed':>9}{'mixed %':>9}")
    total_pairs = total_mixed = 0
    empty_scenes = []
    for scene_path in sorted(p for p in run_dir.iterdir() if p.is_dir()):
        database_path = scene_path / "database.db"
        if not database_path.exists():
            print(f"{scene_path.name:<20}{'-- no database --':>34}")
            continue
        with sqlite3.connect(f"file:{database_path}?mode=ro", uri=True) as db:
            is_calibrated = {
                camera_id: prior != 0
                for camera_id, prior in db.execute(
                    "SELECT camera_id, prior_focal_length FROM cameras"
                )
            }
            camera_of_image = dict(
                db.execute("SELECT image_id, camera_id FROM images")
            )
            pair_ids = [
                pair_id
                for (pair_id,) in db.execute(
                    "SELECT pair_id FROM two_view_geometries WHERE rows > 0"
                )
            ]
        mixed = 0
        for pair_id in pair_ids:
            image_id2 = pair_id % K_MAX_NUM_IMAGES
            image_id1 = (pair_id - image_id2) // K_MAX_NUM_IMAGES
            calibrated1 = is_calibrated.get(camera_of_image.get(image_id1))
            calibrated2 = is_calibrated.get(camera_of_image.get(image_id2))
            if calibrated1 is not None and calibrated1 != calibrated2:
                mixed += 1
        if not pair_ids:
            empty_scenes.append(scene_path.name)
        num_calibrated = sum(is_calibrated.values())
        percent = 100.0 * mixed / len(pair_ids) if pair_ids else 0.0
        total_pairs += len(pair_ids)
        total_mixed += mixed
        print(
            f"{scene_path.name:<20}{len(is_calibrated):>9}{num_calibrated:>8}"
            f"{len(pair_ids):>9}{mixed:>9}{percent:>8.1f}%"
        )
    percent = 100.0 * total_mixed / total_pairs if total_pairs else 0.0
    print(f"{'TOTAL':<20}{'':>9}{'':>8}{total_pairs:>9}"
          f"{total_mixed:>9}{percent:>8.1f}%")
    if total_mixed == 0:
        print(
            "\nNo mixed pairs: the one-sided focal path was never taken "
            "and the two arms\ncannot differ. Check "
            "--calibrated-fraction, and use a fresh --run-name if the\n"
            "workspaces were built by a differently configured run."
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


def parse_split_args(argv: list[str]) -> argparse.Namespace:
    """Parses the shim-only invocation that primes one scene's cameras."""
    parser = argparse.ArgumentParser(
        prog=f"{Path(__file__).name} --set-calibration-split",
        description=set_calibration_split.__doc__,
    )
    parser.add_argument("--set-calibration-split", action="store_true")
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--image-path", type=Path, required=True)
    parser.add_argument("--calibrated-fraction", type=float, required=True)
    parser.add_argument("--split-seed", type=int, required=True)
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
        "--calibrated-fraction",
        type=float,
        default=0.5,
        help="Fraction of each scene's cameras given their GT focal. 0.5 "
        "maximises the number of mixed pairs, which is what the one-sided "
        "focal estimator runs on.",
    )
    parser.add_argument(
        "--split-seed",
        type=int,
        default=0,
        help="Seed for choosing which cameras are calibrated. Independent "
        "of the reconstruction seeds, so both arms and every seed share "
        "one partition. Vary it to check that the result is not specific "
        "to a single split.",
    )
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
        default="one-sided-focal",
        help="Subdirectory of --run-path shared by both arms. Use a fresh "
        "name whenever the configuration changes, since a workspace caches "
        "the database it was built with.",
    )
    parser.add_argument(
        "--scenes",
        nargs="+",
        default=[],
        metavar="SCENE",
        help="Scenes to evaluate. Defaults to every ETH3D DSLR scene: "
        "the calibration split is per image, so a scene with several GT "
        "cameras needs no special handling and nothing is excluded.",
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
    if not 0.0 < args.calibrated_fraction < 1.0:
        parser.error(
            "--calibrated-fraction must lie strictly between 0 and 1; at "
            "0 or 1 every pair has two equally (un)calibrated sides and "
            "the one-sided focal path is never taken"
        )
    return args


def main() -> None:
    if "--set-calibration-split" in sys.argv[1:]:
        split_args = parse_split_args(sys.argv[1:])
        set_calibration_split(
            split_args.database,
            split_args.image_path,
            split_args.calibrated_fraction,
            split_args.split_seed,
        )
        return

    args = parse_args()

    if args.compare_only:
        compare_arms(args)
        return

    for path in (args.baseline_colmap, args.branch_colmap):
        if not path.is_file():
            raise SystemExit(f"Not a colmap binary: {path}")
    if args.baseline_colmap.resolve() == args.branch_colmap.resolve():
        print("WARNING: both arms use the same binary; A-B measures only "
              "run-to-run noise.")

    scenes = " ".join(args.scenes) if args.scenes else f"all ETH3D {CATEGORY}"
    print(f"Scenes: {scenes}")
    print(f"Seeds: 0..{args.num_seeds - 1}")
    print(f"Calibrated fraction: {args.calibrated_fraction} "
          f"(split seed {args.split_seed})")
    print(f"Run directory: {args.run_path / args.run_name}")

    # Sequentially: the arms share the scene workspaces and would corrupt each
    # other's databases if run at the same time.
    run_arm(args, args.baseline_label, args.baseline_colmap, args.scenes)
    run_arm(args, args.branch_label, args.branch_colmap, args.scenes)

    summarize_databases(args)
    report_runtime(
        args.run_path / args.run_name,
        (args.branch_label, args.baseline_label),
    )
    compare_arms(args)


if __name__ == "__main__":
    main()
