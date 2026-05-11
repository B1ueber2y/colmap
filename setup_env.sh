#!/bin/bash
# COLMAP-IMU dependencies: system 2025-06 stack (gcc@12.2.0) + user's own spack
# Note: boost from system stack (lfe3zyd...) — user's spack boost lacks graph/program_options components

S=/cluster/software/stacks/2025-06/linux-ubuntu22.04-x86_64_v3/gcc-12.2.0
A=/cluster/home/ademirtas/spack/2025-06/install_tree/linux-ubuntu22.04-x86_64_v3/gcc-12.2.0
CUDA=/cluster/software/stacks/2025-06/linux-ubuntu22.04-x86_64_v3/gcc-12.2.0/cuda-12.6.2-6ovgc2ircbnyntakw37gaj2j7k2fuwp3
CUDNN=$S/cudnn-9.2.0.82-12-epjrq6jaunushvs5vefkbiy5euxwggv4

export PATH=\
$S/cmake-3.30.5-mmkob5m73rcew67mxbgj3gkpvj24jj2f/bin:\
$S/ninja-1.12.1-bqsotl4ftyf6cdfydsfe7ynkuvnr3vvc/bin:\
$CUDA/bin:\
$PATH

export CMAKE_PREFIX_PATH=\
$S/openblas-0.3.28-qfx6vgjp2xp3iilpmzsj2ykjlua7oeut:\
$S/glew-2.2.0-s3yn2hdgkdpvdvuul7fejqz3wwg5llkn:\
$S/ceres-solver-2.2.0-dw474ukdz26qgi7xv5umqippbgxno3ra:\
$S/glog-0.7.1-dpbezgnzfie32r7kvkebvi4mnqwbmucm:\
$S/gflags-2.2.2-i664jziotkkq22fdv2o7fim2d23fj3oz:\
$S/suite-sparse-7.7.0-z72pwkvsiqt4zb3gshli3tyd7p6babkq:\
$S/boost-1.73.0-lfe3zydkav3q5w6opfedkunb3scuqijh:\
$S/eigen-3.4.0-lh4h7c3wed6mhhlsurszucjzd22gmo4y:\
$S/cgal-6.0.1-opwegd7c43ay3xafuxt4x4mjytwfvzep:\
$S/sqlite-3.46.0-q75jhdf4n642ywj6mwdqix3l6fri2l24:\
$S/curl-8.10.1-t4u6ugpvwnggmvs4sorzbcuwz4zg73f6:\
$S/openssl-3.4.0-nxwzaxjdmghh4uqhdyv7nqn6idztcs4a:\
$S/metis-5.1.0-wxm2hb2u6e4m2sudjyagpqtqa4lbkcox:\
$S/gmp-6.3.0-iichog46e24hp4bfppcv2it26wo65an3:\
$S/mpfr-4.2.1-enygljlwsm57ypuqp6yeg4ailpdorsa4:\
$S/mesa-glu-9.0.2-64i4rzmaopawgynfbgh6xkgo4uctpo73:\
$A/openimageio-2.5.15.0-3xoqg6p5gbmsxswanrdzdkuac2tc7jac:\
$A/qt-5.15.15-4ou4zbamxzvs44vc72j7xjhy55ntnwwd:\
$CUDA:\
${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}

export LD_LIBRARY_PATH=\
$S/openblas-0.3.28-qfx6vgjp2xp3iilpmzsj2ykjlua7oeut/lib:\
$S/glew-2.2.0-s3yn2hdgkdpvdvuul7fejqz3wwg5llkn/lib:\
$S/ceres-solver-2.2.0-dw474ukdz26qgi7xv5umqippbgxno3ra/lib:\
$S/glog-0.7.1-dpbezgnzfie32r7kvkebvi4mnqwbmucm/lib:\
$S/gflags-2.2.2-i664jziotkkq22fdv2o7fim2d23fj3oz/lib:\
$S/metis-5.1.0-wxm2hb2u6e4m2sudjyagpqtqa4lbkcox/lib:\
$S/gmp-6.3.0-iichog46e24hp4bfppcv2it26wo65an3/lib:\
$S/mpfr-4.2.1-enygljlwsm57ypuqp6yeg4ailpdorsa4/lib:\
$S/mesa-glu-9.0.2-64i4rzmaopawgynfbgh6xkgo4uctpo73/lib:\
$S/suite-sparse-7.7.0-z72pwkvsiqt4zb3gshli3tyd7p6babkq/lib:\
$S/boost-1.73.0-lfe3zydkav3q5w6opfedkunb3scuqijh/lib:\
$S/sqlite-3.46.0-q75jhdf4n642ywj6mwdqix3l6fri2l24/lib:\
$S/curl-8.10.1-t4u6ugpvwnggmvs4sorzbcuwz4zg73f6/lib:\
$S/openssl-3.4.0-nxwzaxjdmghh4uqhdyv7nqn6idztcs4a/lib64:\
$A/openimageio-2.5.15.0-3xoqg6p5gbmsxswanrdzdkuac2tc7jac/lib:\
$A/qt-5.15.15-4ou4zbamxzvs44vc72j7xjhy55ntnwwd/lib:\
$S/imath-3.1.11-okm224eptfsgvgy55k4vroiab6cj2ao4/lib:\
$CUDA/lib64:\
$CUDNN/lib:\
${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
