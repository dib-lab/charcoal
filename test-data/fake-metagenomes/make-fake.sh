#! /bin/bash
sourmash compute -k 31 --track-abundance --scaled 2000 --merge 2+47+63 ../genomes/{2,47,63}.fa -o 2+47+63.sig
sourmash compute -k 31 --track-abundance --scaled 2000 --merge 2+63 ../genomes/{2,63}.fa -o 2+63.sig
sourmash compute -k 31 --track-abundance --scaled 2000 --merge 47+63 ../genomes/{47,63}.fa -o 47+63.sig
sourmash compute -k 31 --track-abundance --scaled 2000 --merge 2+47 ../genomes/{2,47}.fa -o 2+47.sig
sourmash compute -k 31 --track-abundance --scaled 2000 --merge 1+3+4 ~/dev/sourmash/podar-ref/{1,3,4}.fa -o 1+3+4.sig
