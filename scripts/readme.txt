You need
1. posgres installed
2. pgloader installed
3. a superuser defined for postgres 
    - In your postgres (type psql postgres into your terminal... should work) type:
    - obviously change the username and userpassword... don't use a password that is shared with other things
      this could easily be leaked onto github if you aren't careful.
    CREATE ROLE username WITH LOGIN SUPERUSER PASSWORD 'userpassword';

Then Change the following
    1. Change kevbot:kevbotpass to your info in the sqlite_to_postgres.ipynb
    2. Change the locations in the eQ_sqlite_to_posgres.load as well
