from distutils.core import setup, Extension

triemodule = Extension('trie', sources=['triemodule.c','trie.c',"ngg.c","dnautil.c","convert.c", "ahocorasick/aho_queue.c","ahocorasick/aho_trie.c","ahocorasick/ahocorasick.c"],extra_compile_args=["-std=c99","-w"])


setup (name = 'trie',
       version = '1.0',
       description = 'CRISPR-project trie',
       ext_modules = [triemodule]
       )
