from setuptools import find_packages, setup

setup(
    name='ykktmm',
    packages=find_packages(include=['ykktmm']),
    version='1.0.0',
    description='TMM pra engenharia acustica',
    install_requires=[
        'pandas==0.23.3',
        'numpy>=1.14.5',
        'matplotlib>=2.2.0'],
    author='Bruno Yukio Miyata Silva',
)