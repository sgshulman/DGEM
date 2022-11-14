import os
import sys
from shutil import copyfile
import subprocess as sp

eps: float = 1e-10


def compute_difference(reference: str, filename: str) -> float:
    process = sp.run(["../../utils/imdiff", "-c", reference, filename], capture_output=True)
    refsum_line = process.stdout.decode("utf-8").split("\n")[1].strip()
    return float(refsum_line.split()[1])


def check_difference(reference: str, filename: str) -> bool:
    difference: float = compute_difference(reference, filename)
    if difference > eps:
        print("difference between " + reference + " and " + filename + " is " + str(difference))
        return False
    return True


def run_single_scattering_test(folder: str, observer_pos: str) -> int:
    print()
    print()
    print("Single scattering test " + folder + " started")
    os.chdir(folder)
    copyfile("../../dgem", "./dgem")
    os.chmod("./dgem", 0o775)
    sp.run(["./dgem"])

    ok: bool = check_difference("fimage" + observer_pos + "_00_ref.dat", "fimage" + observer_pos + "_00.dat")
    ok = check_difference("fimage" + observer_pos + "_01_ref.dat", "fimage" + observer_pos + "_01.dat") and ok
    ok = check_difference("qimage" + observer_pos + "_00_ref.dat", "qimage" + observer_pos + "_00.dat") and ok
    ok = check_difference("qimage" + observer_pos + "_00_ref.dat", "qimage" + observer_pos + "_01.dat") and ok
    ok = check_difference("uimage" + observer_pos + "_00_ref.dat", "uimage" + observer_pos + "_00.dat") and ok
    ok = check_difference("uimage" + observer_pos + "_00_ref.dat", "uimage" + observer_pos + "_01.dat") and ok
    os.chdir("../")

    if not ok:
        print("Single scattering test " + folder + " failed")
        return 1

    print("Single scattering test " + folder + " succeeded")
    return 0


def run_tests() -> int:
    failed_number: int = run_single_scattering_test("ps_ss_disk_0_dgem_h_20_point", "00_00")
    failed_number = failed_number + run_single_scattering_test("ps_ss_disk_0_mc_10_5", "00_00")
    failed_number = failed_number + run_single_scattering_test("ps_ss_disk_45_dgem_5_line", "00_45")
    failed_number = failed_number + run_single_scattering_test("ps_ss_disk_45_sobol_10_5", "00_45")
    failed_number = failed_number + run_single_scattering_test("ps_ss_fractal_dgem_5_10", "00_60")
    failed_number = failed_number + run_single_scattering_test("ps_ss_fractal_dgem_h_20_hex", "00_60")
    failed_number = failed_number + run_single_scattering_test("ps_ss_sphere_dgem_5_point", "00_00")
    return failed_number


if __name__ == "__main__":
    failed_number: int = run_tests()
    if failed_number == 0:
        print("All tests passed")
        sys.exit(0)

    print(str(failed_number) + " tests failed")
    sys.exit(1)
