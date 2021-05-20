# syntax = docker/dockerfile:1.2
FROM continuumio/miniconda3:4.9.2
SHELL ["/bin/bash", "-c"]

WORKDIR /pipeline

COPY covid-pipeline.yml .

RUN --mount=type=cache,target=/opt/conda/pkgs conda env create --name covidpipeline -f covid-pipeline.yml

COPY requirements.txt .

RUN --mount=type=cache,target=/root/.cache source /opt/conda/etc/profile.d/conda.sh && \
    conda activate covidpipeline && \
    pip install -r requirements.txt && \
    pip install multiqc

RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda activate covidpipeline && \
    npm install --global @neherlab/nextclade

RUN apt update && apt install libncurses5 -y

COPY . .

RUN chmod +x pipeline.py

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "covidpipeline", "python", "/pipeline/pipeline.py", "all" ]

CMD [ "--help" ]
