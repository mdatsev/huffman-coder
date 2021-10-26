## Compression speed
```
$ time (./a.out -cf test.huf random-1GB)
real    1m59.529s
user    1m7.191s
sys     0m8.608s
```

## Decompression speed
```
$ time (./a.out -xf test.huf random-1GB)
real    0m31.712s
user    0m15.883s
sys     0m1.995s
```

## Compression space efficiency
```
tests/one-letter-file-big ( 1063335 bytes saved; 88% )
tests/bulgarian ( 6507 bytes saved; 50% )
tests/latin ( 5244 bytes saved; 47% )
tests/english ( 3084 bytes saved; 43% )
tests/code.cpp ( 16871 bytes saved; 41% )
a.out ( 185185 bytes saved; 23% )
tests/random-1MB ( 0 bytes saved; 0% )
```