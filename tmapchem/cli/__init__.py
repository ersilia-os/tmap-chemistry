import click

@click.command()
@click.option('--input_file', '-i', help='Input file in CSV format')
@click.option("--output_folder", '-o', default=None, help='Output_folder')
@click.option('--name', '-n', default=None, help='Name of the visualization')
def cli(input_file, output_folder, name):
	print(input_file, output_folder, name)


if __name__ == "__main__":
    cli()