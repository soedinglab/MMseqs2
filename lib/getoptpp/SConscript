Import ('env')

name = 'getoptpp'
inc = env.Dir('.')
ext_inc = env.Dir('getoptpp')
src = env.Glob('src/*.cpp')
deps = []

env.CreateStaticLibrary(name, inc, ext_inc, src, deps)

