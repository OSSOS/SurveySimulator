# setup.py
from distutils.dir_util import remove_tree
import pathlib
from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py
import shutil
import subprocess
import sys

# Remove the build directory at the start

if pathlib.PosixPath("build").is_dir():
    remove_tree('build', verbose=True)


ROOT = pathlib.Path(__file__).parent.resolve()
FORTRAN_DIR = ROOT / "F95"
make = shutil.which('make')
ext_name = "_ossssimlib"

# setuptools validates package dirs during egg_info, before any build command
# has run make. Seed an empty ossssimlib package; f90wrap regenerates the real
# contents during build_fortran().
_generated_pkg = FORTRAN_DIR / "ossssimlib"
_generated_pkg.mkdir(exist_ok=True)
(_generated_pkg / "__init__.py").touch(exist_ok=True)

_fortran_built = False


def build_fortran():
    """
    Run the F95 Makefile once per setup.py invocation.

    Produces F95/_ossssimlib.<abi>.so and the f90wrap-generated
    F95/ossssimlib/ package. Both build_py (packages) and build_ext
    (extension) need these, and either may run first.
    """
    global _fortran_built
    if _fortran_built:
        return
    # Drop stale objects/mods and f2py/meson tree before every extension build.
    # Stale debug.mod previously left f90wrap wrappers with ambiguous names.
    subprocess.check_call([make, "clean"], cwd=FORTRAN_DIR)
    try:
        subprocess.check_call([make, "MODULE=ossssimlib"], cwd=FORTRAN_DIR)
    except subprocess.CalledProcessError:
        log = FORTRAN_DIR / "f2py_f90wrap.log"
        if log.is_file():
            lines = log.read_text(errors="replace").splitlines()
            print(f"==== {log} (last 100 lines) ====", file=sys.stderr)
            print("\n".join(lines[-100:]), file=sys.stderr)
        raise
    _fortran_built = True


class BuildPyWithMake(build_py):
    def run(self):
        # F95/ossssimlib must exist before packages are collected/copied.
        build_fortran()
        super().run()


class BuildExtWithMake(build_ext):
    def run(self):
        build_fortran()
        super().run()

    def build_extension(self, ext):
        # The Makefile already produced the shared library; just move it to
        # wherever setuptools expects this extension (build-lib for regular
        # installs, the source tree for editable installs).
        assert ext_name == ext.name
        dest = pathlib.Path(self.get_ext_fullpath(ext.name))
        dest.parent.mkdir(parents=True, exist_ok=True)
        src = FORTRAN_DIR / dest.name
        if dest.exists():
            dest.unlink()
        # shutil.move (not Path.rename): pip's build-lib may live on a
        # different filesystem (e.g. /tmp vs network home on CANFAR), where
        # rename fails with EXDEV "Invalid cross-device link".
        shutil.move(str(src), str(dest))


setup(
    cmdclass={"build_ext": BuildExtWithMake, "build_py": BuildPyWithMake},
)
