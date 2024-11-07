import os

from fastapi import FastAPI

app = FastAPI()


@app.get("/")
async def read_root():
    return {"message": os.getenv("NOT_SECRET_VARIABLE", "Value not provided")}
