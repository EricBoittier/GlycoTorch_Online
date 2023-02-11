import sys
sys.path.append("/Users/ericboittier/PycharmProjects/GlycoTorch_Online")
from GlycoTorch import app

print("app.main: __name__ is", __name__)
print(app.__dict__)

# if __name__ == '__main__':
    # app.run(debug=True)