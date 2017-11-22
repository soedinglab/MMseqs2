class Mmseqs2 < Formula
  desc "Software suite for very fast protein sequence search and clustering"
  homepage "https://mmseqs.org"
  # tag "bioinformatics"
  # doi "10.1038/nbt.3988"
  url "https://github.com/soedinglab/MMseqs2.git", :using => :git, :tag => "1-c7a89", :revision => "c7a89b6aa6d94796783e94d49750b3e0882573aa"
  version "1-c7a89"
  head do
    url "https://github.com/soedinglab/MMseqs2.git", :using => :git
    resource "documentation" do
      url "https://github.com/soedinglab/MMseqs2.wiki.git", :using => :git
    end
  end

  option "with-avx2", "Build MMseqs2 with AVX2 instruction set support"
  option "with-sse41", "Build MMseqs2 with SSE4.1 instruction set support"

  needs :openmp
  depends_on :mpi => [:cxx, :optional]

  depends_on "cmake" => :build
  depends_on :arch => :x86_64

  resource "documentation" do
    url "https://github.com/soedinglab/MMseqs2.wiki.git", :using => :git, :revision => "6dbd3666edb64fc71173ee714014e88c1ebe2dfc"
  end

  def install
    args = *std_cmake_args

    args << (build.with?("mpi") ? "-DHAVE_MPI=1" : "-DHAVE_MPI=0")
    args << (build.with?("test") ? "-DHAVE_TESTS=1" : "-DHAVE_TESTS=0")

    if build.bottle?
      args << "-DHAVE_SSE4_1=1"
    else
      args << "-DHAVE_AVX2=1" if build.with?("avx2")
      args << "-DHAVE_SSE4_1=1" if build.with?("sse41")
    end

    system "cmake", ".", *args, "--debug-output", "--trace"
    system "make"

    bin.install "src/mmseqs"
    doc.install "README.md"
    doc.install resource("documentation")
    doc.install "examples"
    bash_completion.install "util/bash-completion.sh" => "mmseqs.sh"
  end

  def caveats
    s = "MMseqs2 requires at least SSE4.1 CPU instruction support:\n"
    if Hardware::CPU.sse4?
      s += "SSE4.1 ✔. MMseqs2 will work correctly.\n"
    else
      s += "SSE4.1 ✘. MMseqs2 will not work.\n"
    end
    s += "\nIf supported, AVX2 CPU instructions can be used:\n"
    if Hardware::CPU.avx2?
      s += "AVX2 ✔. MMseqs2 will use AVX2.\n"
    else
      s += "AVX2 ✘. MMseqs2 will fall back to SSE4.1."
    end
    s
  end

  test do
    system "#{bin}/mmseqs", "createdb", "#{doc}/examples/QUERY.fasta", "q"
    system "#{bin}/mmseqs", "cluster", "q", "res", "tmp", "-s", "1", "--cascaded"
  end
end
