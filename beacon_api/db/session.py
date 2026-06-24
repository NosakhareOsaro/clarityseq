"""
beacon_api.db.session
=====================
Async PostgreSQL session factory for the Beacon v2.1.1 API.

Uses SQLAlchemy 2.x async engine with asyncpg driver.
Connection parameters are loaded from environment variables.

Environment variables:
    BEACON_DB_URL: PostgreSQL connection URL
        (e.g. ``postgresql+asyncpg://user:pass@host:5432/beacon``).
        Default: ``postgresql+asyncpg://beacon:beacon@localhost:5432/beacon``.
    BEACON_DB_POOL_SIZE: Connection pool size (default 5).
    BEACON_DB_MAX_OVERFLOW: Pool max overflow connections (default 10).
    BEACON_DB_POOL_TIMEOUT: Pool connection timeout seconds (default 30).

Why asyncpg:
    asyncpg provides native async PostgreSQL access with excellent performance
    (benchmarks show ~3x faster than psycopg2 for async workloads).
    It integrates with SQLAlchemy 2.x via the ``postgresql+asyncpg`` dialect.

References:
    SQLAlchemy 2.x async docs: https://docs.sqlalchemy.org/en/20/orm/extensions/asyncio.html
    asyncpg: https://magicstack.github.io/asyncpg/
"""

from __future__ import annotations

import logging
import os
from collections.abc import AsyncGenerator
from typing import Any

from sqlalchemy.ext.asyncio import (
    AsyncSession,
    async_sessionmaker,
    create_async_engine,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration from environment variables
# ---------------------------------------------------------------------------

_DB_URL: str = os.getenv(
    "BEACON_DB_URL",
    "postgresql+asyncpg://beacon:beacon@localhost:5432/beacon",
)

_POOL_SIZE: int = int(os.getenv("BEACON_DB_POOL_SIZE", "5"))
_MAX_OVERFLOW: int = int(os.getenv("BEACON_DB_MAX_OVERFLOW", "10"))
_POOL_TIMEOUT: int = int(os.getenv("BEACON_DB_POOL_TIMEOUT", "30"))

# ---------------------------------------------------------------------------
# Async engine and session factory
# ---------------------------------------------------------------------------

# Create the SQLAlchemy 2.x async engine with asyncpg driver.
# pool_pre_ping=True enables connection health checks before each use,
# preventing errors from stale connections in long-running applications.
_engine = create_async_engine(
    _DB_URL,
    pool_size=_POOL_SIZE,
    max_overflow=_MAX_OVERFLOW,
    pool_timeout=_POOL_TIMEOUT,
    pool_pre_ping=True,  # health-check connections before use
    echo=os.getenv("BEACON_DB_ECHO", "false").lower() == "true",
)

# Session factory with expire_on_commit=False to prevent lazy-loading
# errors when session objects are accessed after commit in async context.
_AsyncSessionFactory: async_sessionmaker[AsyncSession] = async_sessionmaker(
    bind=_engine,
    class_=AsyncSession,
    expire_on_commit=False,
    autoflush=False,
    autocommit=False,
)


async def get_session() -> AsyncGenerator[AsyncSession, None]:
    """FastAPI dependency: yield an async SQLAlchemy session.

    Provides a transactional session that is automatically rolled back
    on exception and closed on completion.  Use as a FastAPI ``Depends()``
    injection.

    Yields:
        AsyncSession: An async SQLAlchemy session connected to the Beacon
            PostgreSQL database.

    Raises:
        sqlalchemy.exc.SQLAlchemyError: On database connection or query errors.

    Examples:
        In a FastAPI route::

            from fastapi import Depends
            from beacon_api.db.session import get_session
            from sqlalchemy.ext.asyncio import AsyncSession

            @router.get("/example")
            async def example(session: AsyncSession = Depends(get_session)):
                result = await session.execute(select(BeaconVariant))
                return result.scalars().all()
    """
    async with _AsyncSessionFactory() as session:
        try:
            yield session
            await session.commit()
        except Exception:
            await session.rollback()
            raise
        finally:
            await session.close()


async def create_all_tables() -> None:
    """Create all database tables from the ORM models.

    Runs ``CREATE TABLE IF NOT EXISTS`` for all tables registered with
    the declarative Base.  Called at application startup.

    Returns:
        None.

    Note:
        In production, use Alembic migrations instead of create_all().
        This function is intended for development and testing.
    """
    from beacon_api.db.models import Base  # avoid circular imports

    async with _engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)
    logger.info("Beacon database tables created (or already exist).")


async def dispose_engine() -> None:
    """Dispose the async engine connection pool.

    Called at application shutdown to cleanly close all connections.

    Returns:
        None.
    """
    await _engine.dispose()
    logger.info("Beacon database engine disposed.")


def get_engine() -> Any:
    """Return the global async SQLAlchemy engine.

    Returns:
        AsyncEngine instance.

    Note:
        Prefer using ``get_session()`` for database operations.
        Use this function only for schema management (create_all, Alembic).
    """
    return _engine
