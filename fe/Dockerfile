FROM python:3.12-slim
WORKDIR /app
COPY . /app
RUN pip install --quiet --no-cache-dir -r requirements.txt
EXPOSE 80
ENV PORT=80
CMD ["gunicorn", "-b", "0.0.0.0:80", "app:server"]
