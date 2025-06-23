import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse

# results에서 csv 파일 읽고 그래프 생성하기 위한 코드

def main():
    parser = argparse.ArgumentParser(description = "Plot pore size distribution from CSV")
    parser.add_argument(
        "input_csv",
        help = "path to input CSV file (e.g. results/result_ex1_1.csv)")
    parser.add_argument(
        "output_png",
        help = "path to output PNG file (e.g. results/plots/plot_ex1_1.png)")
    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)

    df.columns = [col.strip() for col in df.columns]

    # total = df['reduced accessible volume'].sum()
    # df['RecuedAccessibleVolume'] = df['reduced accessible volume'] / total
    df['RecuedAccessibleVolume'] = df['reduced accessible volume']

    out_dir = os.path.dirname(args.output_png)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # 그래프
    plt.figure()
    plt.plot(df['diameter(nm)'], df['RecuedAccessibleVolume'], '-o')
    plt.xlabel('Diameter (nm)')
    plt.ylabel('Recued Accessible Volume')
    plt.title('Graphite Pore Size Distribution')
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(args.output_png)

    print(f"Plot saved to {args.output_png}")

if __name__ == "__main__":
    main()