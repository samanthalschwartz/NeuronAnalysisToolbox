 dir = getDirectory("Choose Source Directory ");
 list = getFileList(dir);
 n = list.length;
 print(n)
 j = 1;
 print(list[j])
 // first get extension
 extdelim = "."
 extsplit= split(list[j],extdelim);
 fullname = extsplit[0];
 extension = extsplit[1];
 extension = "." + extension;
 print(extsplit.length)
 print("fullname  = "+fullname)
 print("extension = " + extension)
//find all splits
splitdelims = "(_s|_t|_w)"
splits_split= split(fullname,splitdelims);
print("splits = " + splits_split.length)
print(splits_split[0])
print(splits_split[1])
print(splits_split[2])