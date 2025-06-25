import argparse

def main(name):
    print(name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Say hello')
    parser.add_argument('name', help='your name, enter it')
    args = parser.parse_args()

    main(args.name)

