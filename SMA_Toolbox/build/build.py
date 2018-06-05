#!/usr/bin/env python3
#
# Mark J. Olah [mjo@cs.unm.edu] (03/2014)
#
# Automate builds
#
import os
import os.path
import sys
import argparse
import contextlib


SMA_ROOT=os.path.abspath("..")

BUILD_TOOLS_FILE=os.path.join(SMA_ROOT,'tools','BuildTools.txt')
def get_build_tools():
    with contextlib.closing(open(BUILD_TOOLS_FILE,'r')) as f:
        lines=[l.strip() for l in f.readlines()]
        lines=[l for l in lines if not l.startswith('#')]
        print(lines)
    return lines

BUILD_TOOLS=get_build_tools();
BUILD_TYPES=[ "w64.debug", "w64.release", "linux.debug", "linux.release"]
BUILD_CMAKE={"w64.debug":"W64-Debug.cmake",
             "w64.release":"W64-Release.cmake",
             "linux.debug":"Linux-Debug.cmake",
             "linux.release":"Linux-Release.cmake"}
parser=argparse.ArgumentParser()
parser.add_argument("--jobs","-j", type=int, default=8, help="The number of jobs to run simulataneosly.")
parser.add_argument("--force","-f", action='store_true', help="force the removal and build of directories")
parser.add_argument("--pretend","-p", action='store_true', help="Pretend only.  Do not make any changes.")
parser.add_argument("--verbose","-v", action='store_true', help="Verbose makefile.")
parser.add_argument("--nomake","-n", action='store_true', help="Configure only.  Do not make.")
parser.add_argument("--install","-i", action='store_true', help="Make install")
parser.add_argument("build_type", help="The build type:<%s>"%("|".join(BUILD_TYPES)) )
parser.add_argument("tools",default=None, nargs="*", help="List of tools to build:<%s>"%("|".join(BUILD_TOOLS)) )

RMDIR_CMD="rm -rf %s"
CMAKE_CMD="cd %s && cmake -C ../../../cmake/%s ../../../tools/%s"
MAKE_CMD="make -j8 -C %s JOBS=%i"
MAKE_VERBOSE_CMD="make -j8 -C %s VERBOSE=1 JOBS=%i"
MAKE_INSTALL_CMD="make -C %s install"
MAKE_INSTALL_VERBOSE_CMD="make -C %s VERBOSE=1 install"



def main():
    args=parser.parse_args()
    build_tools=args.tools
    print(build_tools);
    if not build_tools:
        build_tools=BUILD_TOOLS;
    build_type=args.build_type.rstrip(" /")
    build_dir_base=os.path.join(os.getcwd(), build_type)
    build_toolchain=BUILD_CMAKE[build_type]
    print("Build Tools:", build_tools);
    for proj in build_tools:
        build_tool=None;
        for btool in BUILD_TOOLS:
            if proj.lower()==btool.lower():
                build_tool=btool
                break
        if build_tool is None:
            raise("Unable to Build Project: %s",proj)
        build_dir=os.path.join(build_dir_base,build_tool);
        print("Build tool: %s Build type: %s Build toolchain: %s Build dir: %s"%(build_tool, build_type, build_toolchain, build_dir))
        if os.path.exists(build_dir):
            if not args.force:
                print("Build directory: %s Exists."%build_dir)
                ans=input("Force remval of directory (yes/no)?:")
                if ans != "yes":
                    print("Aborting.  No Changes made.")
                    sys.exit(2)
            cmd=RMDIR_CMD%build_dir
            print("Executing: %s"%cmd)
            if not args.pretend:
                os.system(cmd)
        #CMake
        if not args.pretend:
            os.mkdir(build_dir)
        cmd=CMAKE_CMD%(build_dir, build_toolchain, build_tool)
        print("Executing: %s"%cmd)
        if not args.pretend:
            os.system(cmd)
        #Make
        if not args.nomake:
            if args.verbose:
                cmd=MAKE_VERBOSE_CMD%(build_dir, args.jobs)
            else:
                cmd=MAKE_CMD%(build_dir, args.jobs)
            print("Executing: %s"%cmd)
            if not args.pretend:
                os.system(cmd)
            #Make Install
            if args.install:
                if args.verbose:
                    cmd=MAKE_INSTALL_VERBOSE_CMD%(build_dir)
                else:
                    cmd=MAKE_INSTALL_CMD%(build_dir)
                print("Executing: %s"%cmd)
                if not args.pretend:
                    os.system(cmd)



if __name__=="__main__":
    main()
