"""PyInstaller entry point for the TEscope desktop application."""

import multiprocessing
import sys

from tescope.gui.app import main

if __name__ == "__main__":
    multiprocessing.freeze_support()  # required for frozen Windows executables
    sys.exit(main())
