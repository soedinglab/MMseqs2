Summary: FFindex is a very simple index/database for huge amounts of small files. 
Name: ffindex
Version: 0.9.8
Release: 1
License: Create Commons license "Attribution-ShareAlike 3.0"
Group: Utilities/System
Source: http://pubshare.genzentrum.lmu.de/scientific_computing/software/ffindex/ffindex-0.9.8.tar.gz
%description
FFindex is a very simple index/database for huge amounts of small files. The
files are stored concatenated in one big data file, seperated by '\0'. A second
file contains a plain text index, giving name, offset and length of of the
small files. The lookup is currently done with a binary search on an array made
from the index file.

%prep
%setup

%build
make

%install
make install INSTALL_DIR=%{buildroot}%{_prefix}

%files
%doc README LICENSE

/usr/bin/ffindex_build
/usr/bin/ffindex_get
/usr/bin/ffindex_modify
/usr/bin/ffindex_from_fasta
/usr/bin/ffindex_apply
/usr/bin/ffindex_unpack
/usr/include/ffindex.h
/usr/include/ffutil.h
%{_libdir}/libffindex.a
%{_libdir}/libffindex.so
%{_libdir}/libffindex.so.0.1
