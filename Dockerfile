FROM python:3.10.12

COPY . .

RUN pip install wiggelen

RUN pip install gff3-parser
RUN pip install numpy
RUN pip install pyBigWig
RUN pip install pandas
RUN pip install matplotlib
RUN pip install tqdm
RUN pip install scipy