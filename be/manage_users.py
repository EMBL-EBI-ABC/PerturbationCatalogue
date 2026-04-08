#!/usr/bin/env python3
"""CLI tool for managing AI Explorer user accounts.

Usage:
  python manage_users.py create --email user@ebi.ac.uk --name "Display Name"
  python manage_users.py list
  python manage_users.py disable --email user@ebi.ac.uk
  python manage_users.py enable --email user@ebi.ac.uk
  python manage_users.py reset-password --email user@ebi.ac.uk

Reads PG_HOST, PG_PORT, PG_USER, PG_PASSWORD, PG_DB from environment variables
(same as the backend). Set them or use a .env file.
"""

import argparse
import getpass
import os
import sys

import bcrypt
import psycopg2
from dotenv import load_dotenv

load_dotenv()


def get_connection():
    return psycopg2.connect(
        host=os.environ["PG_HOST"],
        port=os.environ["PG_PORT"],
        user=os.environ["PG_USER"],
        password=os.environ["PG_PASSWORD"],
        dbname=os.environ["PG_DB"],
    )


def cmd_create(args):
    password = getpass.getpass("Password: ")
    confirm = getpass.getpass("Confirm password: ")
    if password != confirm:
        print("Passwords do not match.", file=sys.stderr)
        sys.exit(1)
    if len(password) < 8:
        print("Password must be at least 8 characters.", file=sys.stderr)
        sys.exit(1)

    password_hash = bcrypt.hashpw(password.encode(), bcrypt.gensalt()).decode()
    conn = get_connection()
    try:
        with conn.cursor() as cur:
            cur.execute(
                "INSERT INTO users (email, password_hash, display_name) VALUES (%s, %s, %s) RETURNING id",
                (args.email, password_hash, args.name),
            )
            user_id = cur.fetchone()[0]
            conn.commit()
            print(f"Created user {args.email} (id={user_id})")
    except psycopg2.errors.UniqueViolation:
        print(f"User {args.email} already exists.", file=sys.stderr)
        sys.exit(1)
    finally:
        conn.close()


def cmd_list(args):
    conn = get_connection()
    try:
        with conn.cursor() as cur:
            cur.execute(
                "SELECT id, email, display_name, is_active, created_at, last_login FROM users ORDER BY id"
            )
            rows = cur.fetchall()
            if not rows:
                print("No users found.")
                return
            print(f"{'ID':<5} {'Email':<30} {'Name':<20} {'Active':<8} {'Created':<20} {'Last Login':<20}")
            print("-" * 103)
            for row in rows:
                uid, email, name, active, created, last_login = row
                print(
                    f"{uid:<5} {email:<30} {name:<20} {'Yes' if active else 'No':<8} {str(created)[:19]:<20} {str(last_login or '-')[:19]:<20}"
                )
    finally:
        conn.close()


def cmd_disable(args):
    conn = get_connection()
    try:
        with conn.cursor() as cur:
            cur.execute(
                "UPDATE users SET is_active = FALSE WHERE email = %s RETURNING id",
                (args.email,),
            )
            if cur.fetchone():
                conn.commit()
                print(f"Disabled user {args.email}")
            else:
                print(f"User {args.email} not found.", file=sys.stderr)
                sys.exit(1)
    finally:
        conn.close()


def cmd_enable(args):
    conn = get_connection()
    try:
        with conn.cursor() as cur:
            cur.execute(
                "UPDATE users SET is_active = TRUE WHERE email = %s RETURNING id",
                (args.email,),
            )
            if cur.fetchone():
                conn.commit()
                print(f"Enabled user {args.email}")
            else:
                print(f"User {args.email} not found.", file=sys.stderr)
                sys.exit(1)
    finally:
        conn.close()


def cmd_reset_password(args):
    password = getpass.getpass("New password: ")
    confirm = getpass.getpass("Confirm new password: ")
    if password != confirm:
        print("Passwords do not match.", file=sys.stderr)
        sys.exit(1)
    if len(password) < 8:
        print("Password must be at least 8 characters.", file=sys.stderr)
        sys.exit(1)

    password_hash = bcrypt.hashpw(password.encode(), bcrypt.gensalt()).decode()
    conn = get_connection()
    try:
        with conn.cursor() as cur:
            cur.execute(
                "UPDATE users SET password_hash = %s WHERE email = %s RETURNING id",
                (password_hash, args.email),
            )
            if cur.fetchone():
                conn.commit()
                print(f"Password reset for {args.email}")
            else:
                print(f"User {args.email} not found.", file=sys.stderr)
                sys.exit(1)
    finally:
        conn.close()


def main():
    parser = argparse.ArgumentParser(description="Manage AI Explorer user accounts")
    sub = parser.add_subparsers(dest="command", required=True)

    p_create = sub.add_parser("create", help="Create a new user")
    p_create.add_argument("--email", required=True)
    p_create.add_argument("--name", required=True, help="Display name")
    p_create.set_defaults(func=cmd_create)

    p_list = sub.add_parser("list", help="List all users")
    p_list.set_defaults(func=cmd_list)

    p_disable = sub.add_parser("disable", help="Disable a user account")
    p_disable.add_argument("--email", required=True)
    p_disable.set_defaults(func=cmd_disable)

    p_enable = sub.add_parser("enable", help="Enable a user account")
    p_enable.add_argument("--email", required=True)
    p_enable.set_defaults(func=cmd_enable)

    p_reset = sub.add_parser("reset-password", help="Reset a user's password")
    p_reset.add_argument("--email", required=True)
    p_reset.set_defaults(func=cmd_reset_password)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
