
def options(ctx):
    ctx.load("compiler_cxx")

def configure(ctx):
    ctx.load("compiler_cxx")
    # TODO: check provided libraries

def build(ctx):
    ctx.program(
        source=ctx.path.ant_glob('src/**/*.cpp'),
        target='mesher',
        lib='gmsh')