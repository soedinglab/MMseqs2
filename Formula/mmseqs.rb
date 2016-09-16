class Mmseqs < Formula
  desc "MMseqs (Many-against-Many sequence searching) is a software suite for very fast protein sequence searches and clustering of huge protein sequence data sets."
  homepage "https://github.com/soedinglab/mmseqs2"
  head "git@github.com:soedinglab/mmseqs2.git", :using => :git

  option "with-tests", "Build MMseqs test utils"

  depends_on "open-mpi" => :optional
  depends_on "cmake" => :build

  def install
    args = *std_cmake_args

    args << (build.with?("open-mpi") ? "-DHAVE_MPI=1" : "-DHAVE_MPI=0")
    args << (build.with?("tests") ? "-DHAVE_TESTS=1" : "-DHAVE_TESTS=0")

    system "cmake", ".", *args
    system "make"

    bin.install "src/mmseqs"
    doc.install "userguide.pdf"
    bash_completion.install "util/bash-completion.sh" => "mmseqs.sh"
  end

  test do
    system "mmseqs"
  end
end
