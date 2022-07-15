if __name__ == "__main__":
    import sys
    sys.dont_write_bytecode = True
    import descr
    from load import read_task
    task = read_task(sys.argv[1])
    result = descr.run(task)
