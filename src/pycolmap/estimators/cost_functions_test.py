# SPDX-License-Identifier: BSD-3-Clause

import pycolmap


def test_cost_functions_submodule_exists() -> None:
    assert hasattr(pycolmap._core, "cost_functions")


def test_reproj_error_cost_exists() -> None:
    assert hasattr(pycolmap._core.cost_functions, "ReprojErrorCost")


def test_rig_reproj_error_cost_exists() -> None:
    assert hasattr(pycolmap._core.cost_functions, "RigReprojErrorCost")


def test_scaled_rig_reproj_error_cost_exists() -> None:
    assert hasattr(pycolmap._core.cost_functions, "ScaledRigReprojErrorCost")


def test_sampson_error_cost_exists() -> None:
    assert hasattr(pycolmap._core.cost_functions, "SampsonErrorCost")


def test_absolute_pose_prior_cost_exists() -> None:
    assert hasattr(pycolmap._core.cost_functions, "AbsolutePosePriorCost")


def test_absolute_pose_position_prior_cost_exists() -> None:
    assert hasattr(
        pycolmap._core.cost_functions, "AbsolutePosePositionPriorCost"
    )


def test_relative_pose_prior_cost_exists() -> None:
    assert hasattr(pycolmap._core.cost_functions, "RelativePosePriorCost")


def test_point3d_alignment_cost_exists() -> None:
    assert hasattr(pycolmap._core.cost_functions, "Point3DAlignmentCost")


def _make_dummy_preintegrated_data() -> pycolmap.PreintegratedImuData:
    import numpy as np

    calib = pycolmap.ImuCalibration()
    opt = pycolmap.ImuPreintegrationOptions()
    integrator = pycolmap.ImuPreintegrator(opt, calib, 0, 10000000)
    integrator.integrate(
        pycolmap.ImuMeasurement(0, np.zeros(3), np.array([0, 0, 9.81]))
    )
    integrator.integrate(
        pycolmap.ImuMeasurement(10000000, np.zeros(3), np.array([0, 0, 9.81]))
    )
    data = integrator.extract()
    data.finalize()
    return data


def test_visual_centric_imu_preintegration_cost_constructs() -> None:
    cf = pycolmap._core.cost_functions
    data = _make_dummy_preintegrated_data()
    q_iori = pycolmap.Rotation3d()

    c_default = cf.VisualCentricImuPreintegrationCost(data)
    assert c_default is not None

    c_with_q_iori = cf.VisualCentricImuPreintegrationCost(data, q_iori, q_iori)
    assert c_with_q_iori is not None


def test_analytical_visual_centric_imu_preintegration_cost_constructs() -> None:
    cf = pycolmap._core.cost_functions
    data = _make_dummy_preintegrated_data()
    q_iori = pycolmap.Rotation3d()

    c_default = cf.AnalyticalVisualCentricImuPreintegrationCost(data)
    assert c_default is not None

    c_with_q_iori = cf.AnalyticalVisualCentricImuPreintegrationCost(
        data, q_iori, q_iori
    )
    assert c_with_q_iori is not None


def test_inertial_rotation_cost_constructs() -> None:
    cf = pycolmap._core.cost_functions
    data = _make_dummy_preintegrated_data()
    rig = pycolmap.Rigid3d()
    q_iori = pycolmap.Rotation3d()

    c_default = cf.InertialRotationCost(data, rig)
    assert c_default is not None

    c_with_q_iori = cf.InertialRotationCost(data, rig, q_iori, q_iori)
    assert c_with_q_iori is not None


def test_inertial_global_positioning_cost_constructs() -> None:
    cf = pycolmap._core.cost_functions
    data = _make_dummy_preintegrated_data()
    rig = pycolmap.Rigid3d()
    q_cw = pycolmap.Rotation3d()
    q_iori = pycolmap.Rotation3d()

    c_default = cf.InertialGlobalPositioningCost(data, rig, q_cw, q_cw)
    assert c_default is not None

    c_with_q_iori = cf.InertialGlobalPositioningCost(
        data, rig, q_cw, q_cw, q_iori, q_iori
    )
    assert c_with_q_iori is not None


def test_bias_prior_costs_construct() -> None:
    import numpy as np

    cf = pycolmap._core.cost_functions
    prior = np.array([0.01, -0.02, 0.03])

    c_bias = cf.BiasPriorCost(prior, 0.05, 3)
    assert c_bias is not None

    c_gyro = cf.GyroBiasPriorCost(prior, 0.05)
    assert c_gyro is not None

    c_accel = cf.AccelBiasPriorCost(prior, 0.1)
    assert c_accel is not None
