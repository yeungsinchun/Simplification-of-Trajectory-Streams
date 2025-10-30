This is a draft implementation of the paper [Simplification of Trajectory Streams](https://arxiv.org/abs/2503.23025).

To use the code, you should have CMake and CGAL installed on your system. Navigate to build/. Run `cmake ..`, then run `make` to produce the `./simplify` executable. Then run `./simplify < data/taxi_log_2008_by_id/1.txt` or any data of your choosing with the correct format.

Disclaimer: only tested on MacOS (arm64).

/build contains the main executable ./simplify

usage: ./simplify < ../data/taxi_log_2008_by_id

/script
contains the script used to clean the taxi log data to the required format

/data
contains the taxi log data and some artificially created data

Note to myself: fix 6.txt. Result is clearly wrong.

I am grateful that the following paper provides the dataset:
[1] Jing Yuan, Yu Zheng, Xing Xie, and Guangzhong Sun. Driving with knowledge from the physical world.
In The 17th ACM SIGKDD international conference on Knowledge Discovery and Data mining, KDD
’11, New York, NY, USA, 2011. ACM.
[2] Jing Yuan, Yu Zheng, Chengyang Zhang, Wenlei Xie, Xing Xie, Guangzhong Sun, and Yan Huang. Tdrive: driving directions based on taxi trajectories. In Proceedings of the 18th SIGSPATIAL International
Conference on Advances in Geographic Information Systems, GIS ’10, pages 99–108, New York, NY, USA, 2010. ACM.