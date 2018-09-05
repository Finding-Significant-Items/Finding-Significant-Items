# Finding-Significant-Items

## Introduction

Finding top-k frequent items has been a hot issue in data bases. Finding top-k persistent items is a new issue, and has attracted increasing attention in recent years. In practice, users often want to know which items are significant, i.e., not only frequent but also persistent. No prior art can address both of the above two issues at the same time. Also, for high-speed data streams, they cannot achieve high accuracy when the memory is tight. In this paper, we define a new issue, named finding top-k significant items, and propose a novel algorithm namely LTC to handle that issue. LTC can accurately report top-k significant items with tight memory. It includes two key techniques: Long-tail Replacement and a modified CLOCK algorithm. To prove the effectiveness of LTC, we theoretically prove there is no overestimation error and derive the correct rate and error bound. We further conduct extensive experiments on real datasets and synthetic datasets. The experimental results show that LTC achieves 300~10^8 and in average 10^5 times higher accuracy than other related algorithms.

## About the source codes and datasets.

There are three documents. The name of the document means the application of the algorithms in the paper.

For each document, there are several files in it. All mentioned algorithms in the Paper are provided.

The file of LTC.cpp denotes our algorithm, we will detail it below.

We use 3 real datasets and synthetic datasets in our experiments.

Social: This dataset comes from a real social network, which includes users' message and the sending time. We regard the username of the sender of a message as an item ID and the sending time of a message as the timestamp. This dataset contains 1.5M messages, and we divide it into 200 periods with a fixed time interval.

Network: This is a temporal network of interactions on the stack exchange web site. Each item consists of three values u,v,t, which means user u answered user v's question at time t. We regard u as an item ID and t as the timestamp. This dataset contains 10M items, and we divide it into 1000 periods with a fixed time interval. 

CAIDA: This dataset is from CAIDA Anonymized Internet Trace 2016, consisting of IP packets (source IP address, destination IP address, source port, destination port, and protocol type). We regard the source IP address of a packet as an item ID and the index as the timestamp. This dataset contains 10M packets, and we divide it into 500 periods with a fixed time interval.

Synthetic: We generated 5 different datasets according to Zipfian distribution by Web Polygraph, with different skewness from 0.6 to 3.0. We regard the index as the timestamp. Each dataset contains 10M items, and we divide it into 1000 periods at a fixed time interval. We do not conduct experiments on finding persistent items and finding significant items on Synthetic datasets, because every item is distributed randomly, thus the persistent items are always the same as the frequent items on Synthetic datasets.

## How to run

Suppose you've already cloned the repository and installed g++.

Open the terminate and compile the code (e.g., g++ LTC.cpp -o LTC.exe), you will soon get the executable file.

Double-click that executable file, you will be asked to input some parameters (e.g. the number of memory size, k, coefficents and so forth). More details are presented in the file of LTC.cpp.

Here we present a dataset, named stack-new.txt. It contains 10,000,000 flows, and each flow contains the item ID and the arriving time. We commend setting the number of periods to 1000, and the memory size to more than 20KB.

However, since the maximal file size in the Github is 100MB, we have to compress this file. Users should decompress this file if they want to use this dataset.

## Output format

We set some variables in the file of LTC.cpp.

For example, the variable PRE donotes the precision.

There are a lot of notes in the file, and they will give you every information you want.

## Sample Input
10000000
1 1
100
30
1000

where 10000000 means the number of flows, 1 1 means the coefficents of alpha and beta, 100 means the value of k, 30 means the memory size, 1000 means the number of periods.

## Sample Output
precision: 0.79000
ARE: 0.02704