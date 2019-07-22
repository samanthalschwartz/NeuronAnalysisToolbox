var s = "a string";
 macro "Enter String..." {
 s = getString("Enter a String:", s);
 }
 macro "Print String" {
 print("global value of s ("+s+") was set in the first macro");
 } 