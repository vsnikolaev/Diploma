Before running QA fill filelist.
1) At dir with files write

`find $PWD -name "part*.root" | cat > filelist`

2)mv filelist to QA dir

3) run root

`root -l QAstart.c`

