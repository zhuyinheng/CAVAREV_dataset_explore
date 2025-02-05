# CAVAREV_dataset_explore

1. Download images and projection matrix:
    ```sh
    wget https://www5.cs.fau.de/fileadmin/research/datasets/cavarev/cavarev_card.filtered.hq.bin.zip
    unzip cavarev_card.filtered.hq.bin.zip
    wget https://www5.cs.fau.de/fileadmin/research/datasets/cavarev/cavarev.matrices.bin
    wget https://www5.cs.fau.de/fileadmin/research/datasets/cavarev/phasedata_card.txt
    ```

2. Compile and run:
    ```sh
    g++ gated_fdk.cpp cnpy.cpp -o main -lz
    ./main 128 128 128 2 ./cavarev_card.filtered.hq.bin ./cavarev.matrices.bin ./phasedata_card.txt 0.1 0.5 ./reco.vol
    ```

3. The results are saved as `projection_image_0.npy`, etc., which can be visualized using common Python code snippets.

