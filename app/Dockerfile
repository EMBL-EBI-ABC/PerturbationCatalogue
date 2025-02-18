# Use an official Python runtime as a parent image
FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Install any needed packages specified in requirements.txt and gunicorn
RUN pip install --no-cache-dir -r requirements.txt

# Expose port 80
EXPOSE 80

# Define environment variable
ENV PORT=80

# Run the application using gunicorn
CMD ["gunicorn", "-b", "0.0.0.0:80", "app:server"]
