import re
 
foldername = 'Z:\Sam\raw data\190404 Cry2Olig-GephIB-HaloTag unregulated\clone3-1_+BlueLight_timeseries_postblue\'
filename = 'clone3-1_+BlueLight_timeseries_postblue_w1561_s2_t1.TIF'
pattern = re.compile('(.*)_w\d(.*)_s(.*)_t(.*)\.(.*)')
#filename='/Volumes/data/0076-14--2006-01-23/data/--W00088--P00001--Z00000--T00000--nucleus-dapi.tif'
#pattern=re.compile('(.*)--W(.*)--P(.*)--Z(.*)--T(.*)--(.*)\.(.*)')
res = re.search(pattern, filename)
basename = res.group(1)
channel = res.group(2)
position = res.group(3)
time = res.group(4)
print 'prefix', basename
print 'channel', channel
print 'position', position
print 'time', time