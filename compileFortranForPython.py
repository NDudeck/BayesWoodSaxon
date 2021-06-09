from numpy.distutils.core import Extension

ext1 = Extension(name='get_wspot_wf',
                 sources=['get_wspot_schro.f'])

if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(name='f2py_example',
          description="F2PY Users Guide examples",
          author="Pearu Peterson",
          author_email="pearu@cens.ioc.ee",
          ext_modules=[ext1]
          )
