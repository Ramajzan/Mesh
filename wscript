def options(ctx):
    ctx.load('compiler_cxx')

def configure(ctx):
    ctx.load('compiler_cxx')
    ctx.check_cxx(uselib_store='gmsh', fragment='#include "{0}"\nint main() {{ return 0; }}\n'.format('gmsh.h'), lib='gmsh')
    ctx.check_cxx(uselib_store='nglib', fragment='#include <cstddef>\n#include "{0}"\nint main() {{ return 0; }}\n'.format('nglib.h'), lib='nglib')

def build(ctx):
    ctx.program(
        source=ctx.path.ant_glob('src/**/*.cpp'),
        target='mesher',
        use='gmsh nglib tetgen')
