import os
import sys
from colorama import Fore, Style
from shutil import copyfile
import subprocess as sp
import typing

eps: float = 1e-10

single_scattering_tests: typing.List[typing.Tuple[str, str]] = [
    ("ps_ss_disk_0_dgem_h_20_point", "00_00"),
    ("ps_ss_disk_0_mc_10_5", "00_00"),
    ("ps_ss_disk_45_dgem_5_line", "00_45"),
    ("ps_ss_disk_45_sobol_10_5", "00_45"),
    ("ps_ss_fractal_dgem_5_10", "00_60"),
    ("ps_ss_fractal_dgem_h_20_hex", "00_60"),
    ("ps_ss_sphere_dgem_5_point", "00_00")]

multiple_scattering_tests: typing.List[typing.Tuple[str, str]] = [
    ("ps_ms_sphere_mc_10_5", "00_00")]


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
    print("Single scattering test " + Style.BRIGHT + folder + Style.RESET_ALL + " started")
    os.chdir(folder)
    copyfile("../../dgem", "./dgem")
    os.chmod("./dgem", 0o775)
    completed_process : sp.CompletedProcess = sp.run(["./dgem"])
    ok: bool = completed_process.returncode == 0

    if ok:
        ok = check_difference("fimage" + observer_pos + "_00_ref.dat", "fimage" + observer_pos + "_00.dat")
        ok = check_difference("fimage" + observer_pos + "_01_ref.dat", "fimage" + observer_pos + "_01.dat") and ok
        ok = check_difference("qimage" + observer_pos + "_00_ref.dat", "qimage" + observer_pos + "_00.dat") and ok
        ok = check_difference("qimage" + observer_pos + "_00_ref.dat", "qimage" + observer_pos + "_01.dat") and ok
        ok = check_difference("uimage" + observer_pos + "_00_ref.dat", "uimage" + observer_pos + "_00.dat") and ok
        ok = check_difference("uimage" + observer_pos + "_00_ref.dat", "uimage" + observer_pos + "_01.dat") and ok

    os.chdir("../")

    if not ok:
        print(Fore.RED + "Single scattering test " + folder + " failed" + Style.RESET_ALL)
        return 1

    print(Fore.GREEN + "Single scattering test " + folder + " succeeded" + Style.RESET_ALL)
    return 0


def run_multiple_scattering_test(folder: str, observer_pos: str) -> int:
    print()
    print()
    print("Single scattering test " + Style.BRIGHT + folder + Style.RESET_ALL + " started")
    os.chdir(folder)
    copyfile("../../dgem", "./dgem")
    os.chmod("./dgem", 0o775)
    completed_process : sp.CompletedProcess = sp.run(["./dgem"])
    ok: bool = completed_process.returncode == 0

    if ok:
        ok: bool = check_difference("fimage" + observer_pos + "_00_ref.dat", "fimage" + observer_pos + "_00.dat")
        ok = check_difference("fimage" + observer_pos + "_01_ref.dat", "fimage" + observer_pos + "_01.dat") and ok
        ok = check_difference("fimage" + observer_pos + "_02_ref.dat", "fimage" + observer_pos + "_02.dat") and ok
        ok = check_difference("fimage" + observer_pos + "_03_ref.dat", "fimage" + observer_pos + "_03.dat") and ok
        ok = check_difference("fimage" + observer_pos + "_04_ref.dat", "fimage" + observer_pos + "_04.dat") and ok
        ok = check_difference("qimage" + observer_pos + "_01_ref.dat", "qimage" + observer_pos + "_01.dat") and ok
        ok = check_difference("qimage" + observer_pos + "_02_ref.dat", "qimage" + observer_pos + "_02.dat") and ok
        ok = check_difference("qimage" + observer_pos + "_03_ref.dat", "qimage" + observer_pos + "_03.dat") and ok
        ok = check_difference("qimage" + observer_pos + "_04_ref.dat", "qimage" + observer_pos + "_04.dat") and ok
        ok = check_difference("uimage" + observer_pos + "_01_ref.dat", "uimage" + observer_pos + "_01.dat") and ok
        ok = check_difference("uimage" + observer_pos + "_02_ref.dat", "uimage" + observer_pos + "_02.dat") and ok
        ok = check_difference("uimage" + observer_pos + "_03_ref.dat", "uimage" + observer_pos + "_03.dat") and ok
        ok = check_difference("uimage" + observer_pos + "_04_ref.dat", "uimage" + observer_pos + "_04.dat") and ok

    os.chdir("../")

    if not ok:
        print(Fore.RED + "Single scattering test " + folder + " failed" + Style.RESET_ALL)
        return 1

    print(Fore.GREEN + "Single scattering test " + folder + " succeeded" + Style.RESET_ALL)
    return 0


def run_tests() -> typing.Tuple[int, int]:
    failed_ss_tests_number: int = 0

    for test_folder, pos in single_scattering_tests:
        failed_ss_tests_number = failed_ss_tests_number + run_single_scattering_test(test_folder, pos)

    failed_ms_tests_number: int = 0

    for test_folder, pos in multiple_scattering_tests:
        failed_ms_tests_number = failed_ms_tests_number + run_multiple_scattering_test(test_folder, pos)

    return failed_ss_tests_number, failed_ms_tests_number


if __name__ == "__main__":
    print()
    print(Style.BRIGHT + "Start Regression Tests" + Style.RESET_ALL)
    failed_ss_number, failed_ms_number = run_tests()
    print()

    if failed_ss_number == 0 and failed_ms_number == 0:
        print(Style.BRIGHT + Fore.GREEN + "All tests passed!!!" + Style.RESET_ALL)
        sys.exit(0)

    if failed_ss_number == 0:
        print(Style.BRIGHT + Fore.GREEN + "All single scattering tests passed!!!" + Style.RESET_ALL)
    else:
        print(Style.BRIGHT + Fore.RED + str(failed_ss_number) + "/" + str(len(single_scattering_tests))
              + " single scattering tests failed!!!" + Style.RESET_ALL)

    if failed_ms_number == 0:
        print(Style.BRIGHT + Fore.GREEN + "All multiple scattering tests passed!!!" + Style.RESET_ALL)
    else:
        print(Style.BRIGHT + Fore.RED + str(failed_ms_number) + "/" + str(len(multiple_scattering_tests))
              + " multiple scattering tests failed!!!" + Style.RESET_ALL)

    sys.exit(1)
