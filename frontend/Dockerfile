FROM node:20-alpine

WORKDIR /app

COPY frontend/package.json /app/package.json
COPY frontend/tsconfig.json /app/tsconfig.json
COPY frontend/vite.config.ts /app/vite.config.ts
COPY frontend/.eslintrc.cjs /app/.eslintrc.cjs
COPY frontend/index.html /app/index.html
COPY frontend/src /app/src

RUN npm install

EXPOSE 5173

CMD ["npm", "run", "dev", "--", "--host", "0.0.0.0", "--port", "5173"]

