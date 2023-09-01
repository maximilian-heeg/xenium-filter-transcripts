FROM rust AS build
WORKDIR /app
COPY . .
RUN cargo build --release 
RUN strip target/release/xenium-filter-transcripts
ENTRYPOINT ["/app/target/release/xenium-filter-transcripts"]
CMD ["--help"]


FROM debian:bookworm-slim
WORKDIR /app
COPY --from=build /app/target/release/xenium-filter-transcripts /app/
ENTRYPOINT ["/app/xenium-filter-transcripts"]
CMD ["--help"]