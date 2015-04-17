from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "black_scholes",
    ext_modules = cythonize("black_scholes.pyx")
)