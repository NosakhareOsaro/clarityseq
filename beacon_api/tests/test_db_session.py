"""
beacon_api.tests.test_db_session
=================================
pytest tests for the async PostgreSQL session factory in beacon_api.db.session.

Tests cover:
    - get_session() dependency: commit + close on success, rollback + close
      + re-raise on exception.
    - create_all_tables(): delegates to engine.begin() / conn.run_sync().
    - dispose_engine(): delegates to engine.dispose().
    - get_engine(): returns the module-level engine singleton.

The real asyncpg engine is never connected to; SQLAlchemy engine/session
objects are patched with unittest.mock AsyncMock/MagicMock doubles.
"""

from __future__ import annotations

from unittest.mock import AsyncMock, MagicMock, patch

import pytest

from beacon_api.db import session as session_module


class _FakeAsyncSessionCM:
    """Fake async context manager standing in for `_AsyncSessionFactory()`."""

    def __init__(self, session: AsyncMock) -> None:
        self._session = session

    async def __aenter__(self) -> AsyncMock:
        return self._session

    async def __aexit__(self, exc_type, exc, tb) -> bool:
        return False


class _FakeBeginCM:
    """Fake async context manager standing in for `_engine.begin()`."""

    def __init__(self, conn: AsyncMock) -> None:
        self._conn = conn

    async def __aenter__(self) -> AsyncMock:
        return self._conn

    async def __aexit__(self, exc_type, exc, tb) -> bool:
        return False


class TestGetSession:
    """Tests for the get_session() FastAPI dependency generator."""

    @pytest.mark.asyncio
    async def test_success_commits_and_closes(self) -> None:
        """On normal completion, the session is committed then closed."""
        mock_session = AsyncMock()
        fake_factory = MagicMock(return_value=_FakeAsyncSessionCM(mock_session))

        with patch.object(session_module, "_AsyncSessionFactory", fake_factory):
            gen = session_module.get_session()
            yielded = await gen.__anext__()
            assert yielded is mock_session

            with pytest.raises(StopAsyncIteration):
                await gen.__anext__()

        mock_session.commit.assert_awaited_once()
        mock_session.close.assert_awaited_once()
        mock_session.rollback.assert_not_called()

    @pytest.mark.asyncio
    async def test_exception_rolls_back_and_closes(self) -> None:
        """When the caller raises inside the `async with`, rollback + close + re-raise."""
        mock_session = AsyncMock()
        fake_factory = MagicMock(return_value=_FakeAsyncSessionCM(mock_session))

        with patch.object(session_module, "_AsyncSessionFactory", fake_factory):
            gen = session_module.get_session()
            yielded = await gen.__anext__()
            assert yielded is mock_session

            with pytest.raises(ValueError, match="boom"):
                await gen.athrow(ValueError("boom"))

        mock_session.rollback.assert_awaited_once()
        mock_session.close.assert_awaited_once()
        mock_session.commit.assert_not_called()


class TestCreateAllTables:
    """Tests for create_all_tables()."""

    @pytest.mark.asyncio
    async def test_create_all_tables_runs_metadata_create_all(self) -> None:
        """create_all_tables() opens a connection and runs Base.metadata.create_all."""
        mock_conn = AsyncMock()
        fake_begin = MagicMock(return_value=_FakeBeginCM(mock_conn))

        with patch.object(session_module, "_engine") as mock_engine:
            mock_engine.begin = fake_begin
            await session_module.create_all_tables()

        mock_conn.run_sync.assert_awaited_once()


class TestDisposeEngine:
    """Tests for dispose_engine()."""

    @pytest.mark.asyncio
    async def test_dispose_engine_calls_engine_dispose(self) -> None:
        """dispose_engine() awaits _engine.dispose()."""
        with patch.object(session_module, "_engine") as mock_engine:
            mock_engine.dispose = AsyncMock()
            await session_module.dispose_engine()

        mock_engine.dispose.assert_awaited_once()


class TestGetEngine:
    """Tests for get_engine()."""

    def test_get_engine_returns_module_engine(self) -> None:
        """get_engine() returns the module-level _engine singleton."""
        assert session_module.get_engine() is session_module._engine

    def test_get_engine_reflects_patched_engine(self) -> None:
        """get_engine() reflects whatever _engine currently points to."""
        sentinel = object()
        with patch.object(session_module, "_engine", sentinel):
            assert session_module.get_engine() is sentinel


class TestConfiguration:
    """Tests for module-level configuration constants derived from env vars."""

    def test_default_db_url_uses_asyncpg_driver(self) -> None:
        """The default connection URL uses the postgresql+asyncpg dialect."""
        assert "asyncpg" in session_module._DB_URL

    def test_pool_settings_are_ints(self) -> None:
        """Pool size / overflow / timeout are parsed as integers."""
        assert isinstance(session_module._POOL_SIZE, int)
        assert isinstance(session_module._MAX_OVERFLOW, int)
        assert isinstance(session_module._POOL_TIMEOUT, int)
