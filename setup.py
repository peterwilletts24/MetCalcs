from distutils.core import setup

setup(
    name='MetCalcs',
    version='0.1.0',
    author='Peter D. Willetts',
    author_email='peterwilletts24@gmail.com',
    packages=['metcalcs',],
    #scripts=[],    
    license='LICENSE.txt',
    description='Some functions for calculating met parameters'
    long_description=open('README.txt').read()
    install_requires=[
        "numpy",
        "scipy",]
        
