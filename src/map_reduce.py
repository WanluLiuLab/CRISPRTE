# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
# Description: Functions and Methods for functional programming
# Parallizing using multiprocessing.Process

#from pathos.multiprocessing import Pool

from multiprocessing import Process, Manager, Queue
import time

def FLATTEN(x): return [i for s in x for i in s]

class MapReduce(object):
	def __init__(self,  map_func, reduce_func, pool=None, chunk_size=None, num_workers=12):
		super(MapReduce, self).__init__()
		"""
		class MapReduce: implementation of parallel computation using map-reduce
		@args map_func: map function for map-reduce. The function maps an iterable object and returns a mapped list
		@args reduce_func: reduce function for map-reduce. The function reduces the mapped iterable and returns the final result
		@args pool: optional. Possibly Pool from multiprocessing module or pathos.multiprocessing
		@args chunk_size: chunk size for the partition of the list, if provided
		@args num_workers: number of threads to use
		"""
		self.map_func = map_func
		self.reduce_func = reduce_func
		self.chunk_size = chunk_size
		self.Pool = pool

	def partition(self, iterable):
		print("total length: {}".format(len(iterable)))
		return [iterable[x:x+self.chunk_size] for x in range(0, len(iterable), self.chunk_size)]

	def map_func_wrap(self, procnum, return_dict, *args):
		"""
		Wrapper for the map function. A return_dict should be referenced such that the result will be temporally saved
		"""
		return_dict[procnum] = self.map_func(*args)

	def proc_map_reduce(self, iterable, partition=False):
		"""
		Parallel implementation of mapReduce using multiprocessing.Process
		"""
		if partition:
			iterable = self.partition(iterable)
		manager = Manager()
		return_dict = manager.dict()
		proc_queue = []
		for proc,i in enumerate(iterable):
			p = Process(target=self.map_func_wrap, args=(proc, return_dict,i))
			p.Daemon = True
			p.start()
			proc_queue.append(p)
		for p in proc_queue:
			p.join()
		while True:
			if any(proc.is_alive() for proc in proc_queue):
				time.sleep(1)
			else:
				if self.reduce_func:
					return self.reduce_func(return_dict.values())
				else:
					return return_dict.values()

	def proc_map_reduce_test(self, iterable, partition=False):
		"""
		Parallel implementation of mapReduce using multiprocessing.Process
		"""
		if partition:
			iterable = self.partition(iterable)
		manager = Manager()
		return_dict = manager.dict()
		proc_queue = []
		for proc,i in enumerate(iterable):
			f = open("mismatch_proc_{}.txt".format(proc), "w+")
			p = Process(target=self.map_func_wrap, args=(proc, return_dict, f, i))
			p.Daemon = True
			p.start()
			proc_queue.append(p)
		for p in proc_queue:
			p.join(d)
		while True:
			if any(proc.is_alive() for proc in proc_queue):
				time.sleep(1)
			else:
				if self.reduce_func:
					return self.reduce_func(return_dict.values())
				else:
					return return_dict.values()


class ThreadDataList(list):
	def __init__(self, n_threads):
		self.n_threads = n_threads
		self._count = 0
		super(ThreadDataList, self)

	def append(self,*args,**kwargs):
		if self._count >= self.n_threads:
			raise ValueError("Maximum threads exceeds")
		self._count += 1
		super().append(*args,**kwargs)

	def clear(self):
		self._count = 0
		super().clear()