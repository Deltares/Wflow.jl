import argparse
import shutil
import tempfile
import zipfile
from pathlib import Path

import boto3
import tomllib
from botocore.config import Config

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = REPOSITORY_ROOT / "utils/example_models.toml"
INPUT_DIRECTORY = REPOSITORY_ROOT / "Wflow/test/data/input"
BUCKET = "wflow"
ENDPOINT_URL = "https://s3.deltares.nl"
REGION = "eu-west-1"


def add_file(archive: zipfile.ZipFile, source_path: Path, archive_path: str) -> None:
    file_info = zipfile.ZipInfo(archive_path, date_time=(1980, 1, 1, 0, 0, 0))
    file_info.compress_type = zipfile.ZIP_STORED
    file_info.external_attr = 0o644 << 16
    with source_path.open("rb") as source, archive.open(file_info, "w") as target:
        shutil.copyfileobj(source, target)


def create_archive(archive_path: Path, example: dict) -> None:
    config_path = REPOSITORY_ROOT / example["config"]
    if not config_path.is_file():
        raise FileNotFoundError(f"Missing example configuration: {config_path}")

    with zipfile.ZipFile(archive_path, "w") as archive:
        add_file(archive, config_path, config_path.name)
        for input_name in sorted(example["inputs"]):
            input_path = INPUT_DIRECTORY / input_name
            if not input_path.is_file():
                raise FileNotFoundError(f"Missing example input: {input_path}")
            add_file(archive, input_path, f"data/input/{input_name}")


def publish_examples(prefix: str) -> None:
    examples = tomllib.loads(MANIFEST_PATH.read_text(encoding="utf-8"))["examples"]
    client = boto3.client(
        "s3",
        endpoint_url=ENDPOINT_URL,
        region_name=REGION,
        config=Config(s3={"addressing_style": "path"}),
    )

    with tempfile.TemporaryDirectory() as temporary_directory:
        for example_name in sorted(examples):
            example = examples[example_name]
            archive_path = Path(temporary_directory) / example["archive"]
            create_archive(archive_path, example)
            object_name = f"{prefix.strip('/')}/{archive_path.name}"
            print(f"Uploading {archive_path.name} to s3://{BUCKET}/{object_name}")
            client.upload_file(
                str(archive_path),
                BUCKET,
                object_name,
                ExtraArgs={
                    "CacheControl": "public, max-age=300",
                    "ContentType": "application/zip",
                },
            )


def main() -> None:
    parser = argparse.ArgumentParser(description="Publish example model archives")
    parser.add_argument("prefix")
    args = parser.parse_args()
    publish_examples(args.prefix)


if __name__ == "__main__":
    main()
