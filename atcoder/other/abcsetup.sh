#!/bin/bash
#参考リンク：https://iconcreator.hatenablog.com/entry/2021/09/19/200000#fn:3

if [ $# -ne 1 ]; then
  echo "Usage: $0 yxx" 1>&2
  exit 1
fi

if ! [[ "$1" =~ ^[1-9][0-9]{,2}$ ]]; then
  echo "Usage: $0 yxx (x=[0-9],y=[1-9])" 1>&2
  exit 1
fi

if [ ! -e ABC ]; then
  echo "Error: directory 'ABC' not exists." 1>&2
  exit 1
fi

mkdir -p ABC/$1

cp -n template.cpp ABC/$1/a.cpp
cp -n template.cpp ABC/$1/b.cpp
cp -n template.cpp ABC/$1/c.cpp
cp -n template.cpp ABC/$1/d.cpp
cp -n template.cpp ABC/$1/e.cpp
cp -n template.cpp ABC/$1/f.cpp
cp -n template.cpp ABC/$1/g.cpp
cp -n template.cpp ABC/$1/h.cpp

exit 0