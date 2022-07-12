cmake_minimum_required(VERSION 3.11.0)
project(ssim)
find_package(PythonExtensions REQUIRED)

add_library(_surveysub MODULE F95/fortran/surveysub.f95)
python_extension_module(_surveysub)
install(TARGETS _surveysub LIBRARY DESTINATION ssim)