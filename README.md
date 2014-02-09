Targeted L1L2  
====================================

This code implements targeted L1-L2 reconstruction

Installation
------------

The framework requires a recent version of MATLAB.  Some code uses class components that are not expected to work in octave at this time.

Download the code via git

    git clone https://github.com/gregfreeman/targeted_l1l2.git

Download the submodules via git

	cd targeted_l1l2
    git submodule init
    git submodule update

	Download data files  http://research-t-l1l2.s3-website-us-east-1.amazonaws.com/data/refimgs.zip
	unzip images to data/refimgs512x512/ 

Build MEX files

    cd utility/libsvm/matlab
    make
    cd ../../..
    


Getting Started
----------------

Open MATLAB

run the setup script
    setup_targeted_l1l2

run an example test script
    targeted_l1l2_one






