devtools::use_package('rPython')

test_me <- function(){
  rPython::python.load(system.file('python/test.py',package='bcPcaAnalysis'))
  #rPython:::python.call("funcassociate",)
}