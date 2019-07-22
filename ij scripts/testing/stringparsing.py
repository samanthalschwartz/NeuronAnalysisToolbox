
import re
testtxt = 'olig+light-mChcellfill-antiGlyR1000_s5_t1_w3640.TIF';
MJKLab_metamorph_basename = re.compile(r'(_w\d*)|(_s\d*)|(_t\d*)')
mo = MJKLab_metamorph_basename.search(testtxt);
x  = testtxt.split(mo.group(1))
basename = x[0];
MJKLab_metamorph_channel = re.compile(r'(_w\d*)')
MJKLab_metamorph_stagepos = re.compile(r'(_s\d*)')
MJKLab_metamorph_time = re.compile(r'(_t\d*)')
mo = MJKLab_metamorph.findall('olig+light-mChcellfill-antiGlyR1000_s5_t1_w3640.TIF')
print(mo.group(1))

