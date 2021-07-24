# LZJB-2
LZJB-2 compression algorithm
Разработка алгоритма сжатия и декомпрессии для файловой системы. За основу взят алгоритм LZJB. Разработан новый алгоритм, более эффективный.

Описание по ссылке:
https://project-george.blogspot.com/2021/07/lzo-vs-lzjb.html

Сравнение степени сжатия LZ4 vs LZJB-2, размер блока ограничен 4кБайт. Дистанция поиска для алгоритма LZJB-2 =2кБайта. 
Примеры:

сам себя пакует, исполняемый код x86_64

    $ lz4 --best -v -B4096 lzjb2.exe  test.lz4
    using blocks of size 4 KB
    *** LZ4 command line interface 64-bits v1.9.3, by Yann Collet ***
    Compressed 314731 bytes into 157673 bytes ==> 50.10%
    $ lzjb2.exe lzjb2.exe
    ==> 43.14%

небольшой текстовый файл на русском языке в кодировке UTF-8

    $ lz4 --best -v -B4096 test.utf8.txt  test.lz4
    Compressed 5126 bytes into 2895 bytes ==> 56.48%
    $ lzjb2.exe test.utf8.txt
    ==> 47.37%

файл в кодировке base64

    $ lz4 --best -v -B4096 test.b64 test.lz4
    Compressed 466669 bytes into 266322 bytes ==> 57.07%
    $ lzjb2.exe test.b64
    ==> 52.29%

файл векторной графики в формате SVG

    $ lz4 --best -v -B4096 test.svg test.lz4
    Compressed 622684 bytes into 167499 bytes ==> 26.90%
    $ lzjb2.exe test.svg
    ==> 24.59%
