import multiprocessing
from datetime import datetime
from queue import Queue, Empty
import concurrent.futures

def run_parallel_processes(worker_function, data, timeout):
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    queue = Queue()

    def worker(x, return_dict):
        result = worker_function(x[1])
        return_dict[x[0]] = result

    def process_data(i):
        while True:
            try:
                x = queue.get_nowait()
            except Empty:
                break
            else:
                p = multiprocessing.Process(target=worker, args=(x, return_dict))
                p.start()
                p.join(timeout)
                if p.is_alive():
                    print(f"Processing data {x[0]} timed out")
                    p.terminate()
                    p.join()
                else:
                    print(x[0], return_dict[x[0]])
                    continue

    # Put arguments into the queue
    for x in data:
        queue.put(x)

    # Run processes
    with concurrent.futures.ThreadPoolExecutor() as executor:
        executor.map(process_data, range(multiprocessing.cpu_count()))
    
    return return_dict