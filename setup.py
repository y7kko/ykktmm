from setuptools import find_packages, setup

setup(
    name='ykktmm',
    packages=find_packages(include=['ykktmm']),
    version='1.0.0',
    description='TMM pra engenharia acustica',
    install_requires=[
        'numpy>=1.14.0',
        'matplotlib>=2.1.2'],
    author='Bruno Yukio Miyata Silva',
)