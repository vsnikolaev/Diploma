Before  before running flow code fill filelist.
1) cd to directory with root files and write

`find $PWD -name "part*.root" | cat > filelist`

2)mv filelist to Flow directory

3) run root

`root -l Flowstart.c`
