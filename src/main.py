from app import app
from route import *

if __name__ == '__main__':
    app.run(port=5052,  threaded = True)
