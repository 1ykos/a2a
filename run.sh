# g++ -march=native -Wfatal-errors -std=c++17 -ltbb -Og -g -DNDEBUG -I/home/usr/src/partiality/ -I/home/usr/src/asu/ -I/home/usr/src/patchmap/ -I/home/usr/src/wmath/ -I/home/usr/src/geometry/ -I/home/usr/src/dlib/ ./merge.cpp -lboost_container -o merge
./to_easy.awk source <(echo -ne) a2a-overpredict.stream | ./mergeparm2bin > a2a-overpredict.bin
./merge 20 < a2a-overpredict.bin
