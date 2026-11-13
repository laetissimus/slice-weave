"""Entry point of the standalone executable (SliceWeave.exe).

Double-click: starts the local viewer and opens the browser.
From a terminal: the same sub-commands as `sliceweave` (recon, doctor, selftest, simulate, ...).
"""
import multiprocessing
import os
import sys
import traceback

if not getattr(sys, "frozen", False):  # running from the source tree
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _launched_by_double_click() -> bool:
    if os.name != "nt":
        return False
    try:
        import ctypes
        procs = (ctypes.c_uint * 4)()
        # only this process (and the bootloader) attached to the console -> Explorer opened it
        return ctypes.windll.kernel32.GetConsoleProcessList(procs, 4) <= 2
    except Exception:  # noqa: BLE001
        return False


def main():
    multiprocessing.freeze_support()
    pause = _launched_by_double_click()
    args = sys.argv[1:] or ["serve"]
    code = 0
    try:
        from sliceweave.cli import main as cli
        cli(args)
    except SystemExit as e:
        code = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
        if e.code not in (None, 0) and not isinstance(e.code, int):
            print(e.code)
    except Exception:  # noqa: BLE001 - shown to the user before the window closes
        traceback.print_exc()
        code = 1
    if pause and (code != 0 or args[0] != "serve"):
        try:
            input("\nPress Enter to close this window...")
        except EOFError:
            pass
    sys.exit(code)


if __name__ == "__main__":
    main()
